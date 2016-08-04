#!/usr/bin/env python2.7
from __future__ import print_function
import os

from vqsr import gatk_genotype_gvcf
from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.files import upload_or_move_job, get_files_from_filestore

def hard_filter_pipeline(job, uuid, gvcf_id, config):
    """
    Runs GATK Hard Filtering on HaplotypeCaller output

    :param job:
    :param uuid:
    :param vcf_id:
    :param config:
    :return:
    """

    job.fileStore.logToMaster('Running Hard Filter on {}'.format(uuid))
    genotype_gvcf = job.wrapJobFn(gatk_genotype_gvcf, dict(uuid=gvcf_id), config)
    select_snps = job.wrapJobFn(gatk_select_variants, 'SNP',
                                genotype_gvcf.rv(), config)
    snp_filter = job.wrapJobFn(apply_filter, 'SNP', select_snps.rv(), config)
    select_indels = job.wrapJobFn(gatk_select_variants, 'INDEL',
                                  genotype_gvcf.rv(), config)
    indel_filter = job.wrapJobFn(apply_filter, 'INDEL', select_indels.rv(), config)

    job.addChild(genotype_gvcf)
    genotype_gvcf.addChild(select_snps)
    genotype_gvcf.addChild(select_indels)
    select_snps.addChild(snp_filter)
    select_indels.addChild(indel_filter)

    output_dir = os.path.join(config['output_dir'], uuid)
    raw_filename = '%s_raw_variants%s.vcf' % (uuid, config['suffix'])
    snp_filename = '%s_filtered_snps%s.vcf' % (uuid, config['suffix'])
    indel_filename = '%s_filtered_indels%s.vcf' % (uuid, config['suffix'])
    output_raw = job.wrapJobFn(upload_or_move_job, raw_filename, gvcf_id, output_dir=output_dir,
                               ssec=config['ssec'])
    job.addChild(output_raw)
    output_snps = job.wrapJobFn(upload_or_move_job, snp_filename, snp_filter.rv(),
                               output_dir, ssec=config['ssec'])
    output_indels = job.wrapJobFn(upload_or_move_job, indel_filename, indel_filter.rv(),
                                  output_dir, ssec=config['ssec'])
    snp_filter.addChild(output_snps)
    indel_filter.addChild(output_indels)


def gatk_select_variants(job, mode, vcf_id, config):
    """

    :param job:
    :param mode:
    :param vcf_id:
    :param config:
    :return:
    """
    job.fileStore.logToMaster('Running GATK SelectVariants: %s' % mode)
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs['input.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    output = "raw_%s.vcf" % mode

    command = ['-R', 'genome.fa',
               '-T', 'SelectVariants',
               '-V', 'input.vcf',
               '-o', output,
               '-selectType', mode]

    outputs = {output: None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, output))


def apply_filter(job, mode, vcf_id, config):
    """

    :param job:
    :param mode:
    :param vcf_id:
    :param config:
    :return:
    """
    job.fileStore.logToMaster('Apply %s Filter' % mode)
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs['input.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    if mode.upper() == 'SNP':
        expression = '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || ' \
                     'MQRankSum < -12.5 || ReadPosRankSum < -8.0"'

    elif mode.upper() == 'INDEL':
        expression = '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'

    else:
        raise ValueError("Need to specify either SNP or INDEL mode")

    command = ['-T', 'VariantFiltration',
               '-R', 'genome.fa',
               '-V', 'input.vcf',
               '--filterExpression', expression,
               '--filterName', 'GATK_Germline_Hard_Filter_%s' % mode,
               '-o', 'filtered_variants.vcf']

    outputs = {'filtered_variants.vcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'filtered_variants.vcf'))
