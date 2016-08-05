#!/usr/bin/env python2.7
from __future__ import print_function
import os

from toil_scripts.lib.files import get_files_from_filestore
from toil_scripts.lib.programs import docker_call



def run_oncotator(job, config):
    """

    :param job:
    :param bam_id:
    :param bai_id:
    :param vcf_id:
    :param annotations:
    :param config:
    :return:
    """
    job.fileStore.logToMaster('Running Oncotator: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    get_files_from_filestore(job, work_dir, inputs)

    # Call GATK -- HaplotypeCaller
    command = ['-R', 'genome.fa',
               '-T', 'VariantAnnotator',
               '-I', 'input.bam',
               '-V', 'input.vcf',
               '-o', 'output.vcf']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'output.vcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'jpfeil/oncotator:1.9--8fffc356981862d50cfacd711b753700b886b605',
                inputs=inputs.keys(),
                outputs=outputs, mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.vcf'))


def gatk_variant_annotator(job, bam_id, bai_id, vcf_id, annotations, config):
    """

    :param job:
    :param bam_id:
    :param bai_id:
    :param vcf_id:
    :param annotations:
    :param config:
    :return:
    """
    job.fileStore.logToMaster('Running GATK VariantAnnotator: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs['input.bam'] = bam_id
    inputs['input.bam.bai'] = bai_id
    inputs['input.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    # Call GATK -- HaplotypeCaller
    command = ['-R', 'genome.fa',
               '-T', 'VariantAnnotator',
               '-I', 'input.bam',
               '-V', 'input.vcf',
               '-o', 'output.vcf']

    for annotation in annotations:
        command.extend(['-A', annotation])

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'output.vcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs, mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'output.vcf'))
