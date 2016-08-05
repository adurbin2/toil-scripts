#!/usr/bin/env python2.7
from __future__ import print_function
import argparse
import multiprocessing
import os
from copy import deepcopy

import synapseclient
from bd2k.util.processes import which
from bd2k.util.humanize import human2bytes
from toil.job import Job
from urlparse import urlparse
import yaml

from toil_scripts.gatk_germline.germline_config import generate_config, generate_manifest
from toil_scripts.gatk_germline.hard_filter import hard_filter_pipeline
from toil_scripts.gatk_germline.vqsr import vqsr_pipeline
from toil_scripts.lib import require
from toil_scripts.lib.files import get_files_from_filestore, upload_or_move_job
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.urls import download_url_job
from toil_scripts.rnaseq_cgl.rnaseq_cgl_pipeline import generate_file
from toil_scripts.tools.aligners import run_bwakit
from toil_scripts.tools.preprocessing import run_germline_preprocessing, \
    run_samtools_faidx, run_samtools_index, run_samtools_sort, picard_sort_sam, \
    run_picard_create_sequence_dictionary


def run_pipeline(job, samples, config):
    """
    Configures and runs GATK Germline pipeline
    
    :param samples list: List of tuples (uuid, url, rg_line) for joint variant calling.
                         (i.e. samples are grouped together for VQSR). Otherwise,
                         samples are filtered using GATK recommended hard filters.
    :param config dict: Input parameters and reference files
    :return: None

    VQSR and filtering is done if the run_vqsr parameter is set to True. Otherwise, each sample
    will be hard filtered using recommended GATK parateters.
    """
    get_shared_files = Job.wrapJobFn(download_shared_files, config).encapsulate()
    # Pass config dictionary to reference preprocessing
    get_references = Job.wrapJobFn(reference_preprocessing, get_shared_files.rv(),
                                   mock=config['mock_mode']).encapsulate()

    # Run after shared files have been downloaded
    job.addChild(get_shared_files)
    get_shared_files.addChild(get_references)
    # Node to batch each sample job
    batch = Job()
    get_references.addChild(batch)
    # Generate per sample gvcfs {uuid: gvcf_id}
    gvcfs = {}
    for uuid, url, rg_line in samples:
        # get_references.rv() returns an updated version of config
        gvcfs[uuid] = batch.addChildJobFn(gatk_germline_pipeline, uuid, url,
                                           get_references.rv(), rg_line=rg_line).rv()

    if config['run_vqsr'] or len(gvcfs) > 30:
        vqsr_vcf = batch.addFollowOnJobFn(vqsr_pipeline, gvcfs, get_references.rv())

    else:
        for uuid, gvcf in gvcfs.iteritems():
            batch.addFollowOnJobFn(hard_filter_pipeline, uuid, gvcf, get_references.rv())


def gatk_germline_pipeline(job, uuid, url, config, rg_line=None):
    """
    Runs GATK Germline Pipeline on a single sample. Writes gvcf to output directory
    defined in the config dictionary.

    0: Align FASTQ or Download Bam
    1: Generate BAI
    2: GATK Preprocessing (Optional)
    3: GATK HaplotypeCaller
    4: Write GVCF to output directory

    :param job: Toil Job instance
    :param uuid str: Unique identifier for the sample
    :param url str: URL to sample bam or fq file
    :param config dict: Configuration options for pipeline
    :param rg_line str|None: RG line for BWA alignment. Necessary for downstream steps.
    :rtype str: GVCF FileStoreID
    """
    config = deepcopy(config)
    config['uuid'] = uuid
    config['rg_line'] = rg_line

    cores = multiprocessing.cpu_count()
    # Determine how much RAM is available
    if config ['xmx'] is None:
        config['xmx'] = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
    config['xmx'] = human2bytes(config['xmx'])

    # Determine extension on input file
    base, ext1 = os.path.splitext(url)
    _, ext2 = os.path.splitext(base)

    # Produce or download BAM file
    if config['run_bwa'] and rg_line and {ext1, ext2} & {'.fq', '.fastq'}:
        get_bam = job.wrapJobFn(setup_and_run_bwa_kit, url, config).encapsulate()

    elif {ext1, ext2} & {'.bam'}:
        get_bam = job.wrapJobFn(download_url_job, url, name='toil.bam',
                                s3_key_path=config['ssec']).encapsulate()
    else:
        raise ValueError('Accepted formats: FASTQ or BAM\n'
                         'Could not determine {} file type from URL:\n{}'.format(uuid, url))

    sort_bam = job.wrapJobFn(picard_sort_sam, get_bam.rv(), xmx=config['xmx'])
    get_bai = job.wrapJobFn(run_samtools_index, sort_bam.rv())

    job.addChild(get_bam)
    get_bam.addChild(sort_bam)
    sort_bam.addChild(get_bai)

    if config['preprocess']:
        preprocess = job.wrapJobFn(run_germline_preprocessing, cores, uuid,
                                   get_bam.rv(), get_bai.rv(),
                                   config['genome.fa'], config['genome.dict'],
                                   config['genome.fa.fai'], config['phase'],
                                   config['mills'], config['dbsnp'],
                                   mem=config['xmx'], mock=config['mock_mode']).encapsulate()
        get_bai.addChild(preprocess)
        haplotype_caller = job.wrapJobFn(gatk_haplotype_caller, preprocess.rv(0), preprocess.rv(1),
                                         config)

    else:
        haplotype_caller = job.wrapJobFn(gatk_haplotype_caller, get_bam.rv(), get_bai.rv(), config)

    get_bai.addFollowOn(haplotype_caller)
    output_dir = os.path.join(config['output_dir'], uuid)
    vqsr_name = '{}.raw{}.gvcf'.format(uuid, config['suffix'])
    output_gvcf = job.wrapJobFn(upload_or_move_job, vqsr_name, haplotype_caller.rv(),
                                output_dir)
    haplotype_caller.addChild(output_gvcf)
    return haplotype_caller.rv()


def download_shared_files(job, config):
    """
    Downloads reference files shared by all samples in the pipeline

    :param config dict: pipeline configuration options
    :return dict: updated config dictionary with shared fileStoreIDS
    """
    job.fileStore.logToMaster('Downloading shared reference files')
    reference_names = ['genome.fa', 'phase', 'mills', 'dbsnp', 'hapmap', 'omni']
    if config['run_bwa']:
        bwa_references = ['amb', 'ann', 'bwt', 'pac', 'sa', 'alt']
        reference_names.extend(bwa_references)
    for name in reference_names:
        if config[name] == None:
            continue
        try:
            config[name] = job.addChildJobFn(download_url_job, config[name], name=name,
                                             syn=config['syn'], s3_key_path=config['ssec']).rv()
        except KeyError:
            config[name] = None

    return config


def reference_preprocessing(job, config, mock=False):
    """
    Creates a genome.fa.fai and genome.dict file for reference genome.fa
    if one is not already present in the config dictionary

    :param config dict: pipeline configuration options and shared files
                        requires a key called 'genome.fa'
    :return dict: updated config dictionary with reference index files
    """
    job.fileStore.logToMaster('Preparing Reference Files')
    genome_id = config['genome.fa']
    if 'genome.fa.fai' not in config or config['genome.fa.fai'] is None:
        config['genome.fa.fai'] = job.addChildJobFn(run_samtools_faidx, genome_id, mock=mock).rv()
    if 'genome.dict' not in config or config['genome.dict'] is None:
        config['genome.dict'] = job.addChildJobFn(run_picard_create_sequence_dictionary,
                                                  genome_id, mock=mock).rv()
    return config


def setup_and_run_bwa_kit(job, url, config):
    """

    :param job: Toil Job instance
    :param url str: FASTQ sample URL
    :param config dict: Dictionary of configuration parameters
    :return: BAM FileStoreID
    :rtype: str
    """
    job.fileStore.logToMaster("Aligning sample: {}".format(config['uuid']))

    # Download fastq files
    fq1 = job.addChildJobFn(download_url_job, url, name='toil.1.fq', syn=config['syn'],
                            s3_key_path=config['ssec'])
    # Assumes second fastq url is identical
    fq2_url = url.replace('1.fq', '2.fq')
    fq2 = job.addChildJobFn(download_url_job, fq2_url, name='toil.2.fq', syn=config['syn'],
                            s3_key_path=config['ssec'])

    # run_bwakit requires a Namespace object
    bwa_config = argparse.Namespace()
    bwa_config.r1 = fq1.rv()
    bwa_config.r2 = fq2.rv()
    bwa_config.ref = config['genome.fa']
    bwa_config.fai = config['genome.fa.fai']
    keys = ['uuid', 'amb', 'ann', 'bwt', 'pac', 'sa', 'alt', 'rg_line']
    for key in keys:
        setattr(bwa_config, key, config[key])

    cores = multiprocessing.cpu_count()
    return job.addFollowOnJobFn(run_bwakit, bwa_config, cores, sort=True,
                                trim=config['trim']).rv()


def gatk_haplotype_caller(job, bam_id, bai_id, config, annotations=None):
    """
    Uses GATK HaplotypeCaller to identify SNPs and Indels and writes a gVCF.

    :param job: Job instance
    :param bam_id str: bam file store ID
    :param bai_id str: bai file store ID
    :param config dict: pipeline configuration options and shared files
    :return str: gvcf file store ID
    """
    job.fileStore.logToMaster('Running GATK HaplotypeCaller: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs['input.bam'] = bam_id
    inputs['input.bam.bai'] = bai_id
    get_files_from_filestore(job, work_dir, inputs)

    cores = multiprocessing.cpu_count()
    # Call GATK -- HaplotypeCaller
    command = ['-nct', str(cores),
               '-R', 'genome.fa',
               '-T', 'HaplotypeCaller',
               '-I', 'input.bam',
               '-o', 'output.gvcf',
               '-stand_call_conf', '2',
               '-stand_emit_conf', '2',
               '-variant_index_type', 'LINEAR',
               '-variant_index_parameter', '128000',
               '--genotyping_mode', 'Discovery',
               '--emitRefConfidence', 'GVCF']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    if annotations is None:
        annotations = ['QualByDepth', 'FisherStrand', 'StrandOddsRatio',
                       'ReadPosRankSumTest', 'MappingQualityRankSumTest',
                       'RMSMappingQuality']

    for annotation in annotations:
        command.extend(['-A', annotation])

    outputs={'output.gvcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs, mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.gvcf'))


def parse_manifest(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if not line.isspace() and not line.startswith('#'):
                sample = line.strip().split()
                if len(sample) == 2:
                    uuid, url = sample
                    rg_line = None
                    if not url.endswith('bam'):
                        raise ValueError("Expected a BAM file:\n{}:\t{}".format(uuid, url))
                elif len(sample ) == 3:
                    uuid, url, rg_line = sample
                    if not url.endswith('fq') or url.endswith('fastq'):
                        raise ValueError("Expected a FASTA file:\n{}:\t{}".format(uuid, url))
                else:
                    msg = 'Bad manifest format!\n{}'.format(sample)
                    raise ValueError(msg)
                require(urlparse(url).scheme and urlparse(url),
                        'Invalid URL passed for {}'.format(url))
                samples.append((uuid, url, rg_line))
    return samples


def main():
    """
    GATK HaplotypeCaller with variant recalibration.

    Steps in Pipeline:
        1: Download Shared Files
        2: Prepare References
        3: Download Sample
        4: Align FASTQ (Optional)
        5: Call SNPs & INDELs
        6: Genotype and Annotate Variants
        7: Filter Variants

    ===================================================================
    :Dependencies:
    curl            - apt-get install curl
    docker          - apt-get install docker (or 'docker.io' for linux)
    toil            - pip install --pre toil
    """
    # Define Parser object and add to jobTree
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    # Generate subparsers
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the GATK germline pipeline')
    parser_run.add_argument('--config', default='config-toil-germline.yaml', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config".')
    parser_run.add_argument('--manifest', default='manifest-toil-germline.tsv', type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s".')
    parser_run.add_argument('--sample', default=None, nargs=2, type=str,
                       help='Space delimited sample UUID and BAM file in the format: uuid url')
    parser_run.add_argument('--output-dir', default=None, help='Full path to directory or filename where '
                                                               'final results will be output')
    parser_run.add_argument('-s', '--suffix', default=None,
                            help='Additional suffix to add to the names of the output files')
    Job.Runner.addToilOptions(parser_run)
    options = parser.parse_args()

    cwd = os.getcwd()
    if options.command == 'generate-config' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'config-toil-germline.yaml'), generate_config)
    if options.command == 'generate-manifest' or options.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-germline.tsv'), generate_manifest)
    # Pipeline execution
    elif options.command == 'run':
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)), program + ' must be installed on every node.'.format(program))

        require(os.path.exists(options.config), '{} not found. Please run '
                                                '"generate-config"'.format(options.config))
        samples = []
        if options.manifest:
            samples.extend(parse_manifest(options.manifest))
        if options.sample:
            uuid, url, rg_line = options.sample.split()
            samples.append((uuid, url, rg_line))

        # Parse config
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(options.config).read()).iteritems()}
        config_fields = set(config)
        required_fields = {'genome.fa', 'mills', 'dbsnp', 'phase', 'hapmap', 'omni', 'output_dir',
                           'run_bwa', 'preprocess', 'run_vqsr'}
        require(config_fields > required_fields,
                'Missing config parameters:\n{}'.format(', '.join(required_fields - config_fields)))

        if config['output_dir'] is None:
            config['output_dir'] = options.output_dir \
                if options.output_dir else os.path.join(os.getcwd(), 'output')

        if config['suffix'] is None:
            config['suffix'] = options.suffix if options.suffix else ''

        if config['synapse_id'] and config['synapse_pwd']:
            config['syn'] = synapseclient.login(config['synapse_id'], config['synapse_pwd'])
        else:
            config['syn'] = None

        root = Job.wrapJobFn(run_pipeline, samples, config)
        Job.Runner.startToil(root, options)

if __name__ == '__main__':
    main()
