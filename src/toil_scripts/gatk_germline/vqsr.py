#!/usr/bin/env python2.7
from __future__ import print_function
import os
import multiprocessing

from toil_scripts.lib.programs import docker_call
from toil_scripts.lib.files import upload_or_move_job, get_files_from_filestore


def vqsr_pipeline(job, gvcfs, config):
    """
    Takes a dictionary of gvcfs and performs VQSR. Returns filtered vcf

    :param gvcfs:
    :param config:
    :return:
    """
    genotype_gvcf = job.wrapJobFn(gatk_genotype_gvcf, gvcfs, config)
    snp_recal = job.wrapJobFn(gatk_variant_recalibrator_snp, genotype_gvcf.rv(), config)
    indel_recal = job.wrapJobFn(gatk_variant_recalibrator_indel, genotype_gvcf.rv(), config)
    apply_snp_recal = job.wrapJobFn(gatk_apply_variant_recalibration_snp, genotype_gvcf.rv(),
                                    snp_recal.rv(0), snp_recal.rv(1), config)
    apply_indel_recal = job.wrapJobFn(gatk_apply_variant_recalibration_indel, apply_snp_recal.rv(),
                                      indel_recal.rv(), indel_recal.rv(), config)

    job.addFollowOn(genotype_gvcf)
    genotype_gvcf.addChild(snp_recal)
    genotype_gvcf.addChild(indel_recal)
    job.addFollowOn(apply_snp_recal)
    apply_snp_recal.addFollowOn(apply_indel_recal)
    output_dir = os.path.join(config['output_dir'])
    vqsr_name = 'joint.vqsr{}.vcf'.format(config['suffix'])
    output_vqsr = job.wrapJobFn(upload_or_move_job, vqsr_name, apply_indel_recal.rv(),
                                output_dir)
    apply_indel_recal.addChild(output_vqsr)
    return apply_indel_recal.rv()


def gatk_genotype_gvcf(job, gvcf_ids, config, annotations=None):
    """
    Genotypes the gVCF generated by HaplotypeCaller.

    :param job: Job instance
    :param gvcf_ids dict: Dictionary of gvcf file store IDs
    :param config dict: pipeline configuration options and shared files
    :return str: vcf file store ID
    """
    job.fileStore.logToMaster('Running GATK GenotypeGVCF')
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {key: config[key] for key in references}
    inputs.update(gvcf_ids)
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-nt', str(cores),
               '-R', 'genome.fa',
               '-T', 'GenotypeGVCFs',
               '--out', 'joint.vcf',
               '-stand_emit_conf', '10.0',
               '-stand_call_conf', '30.0']

    for uuid in gvcf_ids.keys():
        command.extend(['--variant', uuid])

    if annotations is None:
        annotations = ['QualByDepth', 'FisherStrand', 'StrandOddsRatio', 'ReadPosRankSumTest',
                       'MappingQualityRankSumTest', 'RMSMappingQuality']

    for annotation in annotations:
        command.extend(['-A', annotation])

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs = {'joint.vcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'joint.vcf'))



def gatk_variant_recalibrator_snp(job, vcf_id, config):
    """
    Variant quality score recalibration for SNP variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param config dict: pipeline configuration options and shared files
    :return tuple: recalibration table, tranches file, plots file
    """
    job.fileStore.logToMaster('Running GATK VariantRecalibrator (SNP Mode): {}'.format(config['uuid']))
    cores = multiprocessing.cpu_count()
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict', 'hapmap.vcf', 'omni.vcf', 'dbsnp.vcf', 'phase.vcf']
    inputs = {key: config[key] for key in references}
    inputs['toil.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    # DP is a recommended annotation, but does not work well with exome data
    command = ['-T', 'VariantRecalibrator',
               '-R', 'genome.fa',
               '-input', 'toil.vcf',
               '-nt', str(cores),
               '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', 'hapmap.vcf',
               '-resource:omni,known=false,training=true,truth=false,prior=12.0', 'omni.vcf',
               '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', 'dbsnp.vcf',
               '-resource:1000G,known=false,training=true,truth=false,prior=10.0', 'phase.vcf',
               '-an', 'QD',
               '-an', 'MQ',
               '-an', 'MQRankSum',
               '-an', 'ReadPosRankSum',
               '-an', 'FS',
               '-an', 'SOR',
               '-mode', 'SNP', '-minNumBad', '1000',
               '-recalFile', 'HAPSNP.recal',
               '-tranchesFile', 'HAPSNP.tranches',
               '-rscriptFile', 'HAPSNP.plots']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs = {'HAPSNP.recal': None, 'HAPSNP.tranches': None, 'HAPSNP.plots': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool ='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])

    recal_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPSNP.recal'))
    tranches_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPSNP.tranches'))
    plots_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPSNP.plots'))
    return recal_id, tranches_id, plots_id


def gatk_apply_variant_recalibration_snp(job, vcf_id, recal_id, tranches_id, config):
    """
    Apply variant quality score recalibration for SNP variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param recal_id: recalibration table file store ID
    :param tranches_id: tranches file store ID
    :param config dict: pipeline configuration options and shared files
    :return str: SNP recalibrated VCF file store ID
    """
    job.fileStore.logToMaster('Running GATK ApplyRecalibration (SNP Mode): {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {'toil.vcf': vcf_id,
              'HAPSNP.tranches': tranches_id,
              'HAPSNP.recal': recal_id}
    inputs.update({key: config[key] for key in references})
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-T', 'ApplyRecalibration',
               '-input', 'toil.vcf',
               '-o', 'toil.vqsr.vcf',
               '-R', 'genome.fa',
               '-nt', '1',
               '-ts_filter_level', '99.0',
               '-tranchesFile', 'HAPSNP.tranches',
               '-recalFile', 'HAPSNP.recal',
               '-mode', 'SNP']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'toil.vqsr.vcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'toil.vqsr.vcf'))


def gatk_variant_recalibrator_indel(job, vcf_id, config):
    """
    Variant quality score recalibration for INDEL variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param config dict: pipeline configuration options and shared files
    :return tuple: recalibration table, tranches file, plots file
    """
    job.fileStore.logToMaster('Running GATK VariantRecalibrator (INDEL Mode)')
    cores = multiprocessing.cpu_count()
    work_dir = job.fileStore.getLocalTempDir()
    references = ['genome.fa', 'genome.fa.fai', 'genome.dict', 'mills.vcf']
    inputs = {key: config[key] for key in references}
    inputs['toil.vcf'] = vcf_id
    get_files_from_filestore(job, work_dir, inputs)

    # DP is a recommended annotation, but does not work well with exome data
    command = ['-T', 'VariantRecalibrator',
               '-R', 'genome.fa',
               '-input', 'toil.vcf',
               '-nt', str(cores),
               '-resource:mills,known=true,training=true,truth=true,prior=12.0', 'mills.vcf',
               '-an', 'QD',
               '-an', 'FS',
               '-an', 'SOR',
               '-an', 'ReadPosRankSum',
               '-an', 'MQRankSum',
               '-mode', 'INDEL',
               '-recalFile', 'HAPINDEL.recal',
               '-tranchesFile', 'HAPINDEL.tranches',
               '-rscriptFile', 'HAPINDEL.plots',
               '--maxGaussians', '4',
               '--minNumBadVariants', '5000']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs = {'HAPINDEL.recal': None, 'HAPINDEL.tranches': None, 'HAPINDEL.plots': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool ='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs=outputs,
                mock=config['mock_mode'])
    recal_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPINDEL.recal'))
    tranches_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPINDEL.tranches'))
    plots_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, 'HAPINDEL.plots'))
    return recal_id, tranches_id, plots_id


def gatk_apply_variant_recalibration_indel(job, vcf_id, recal_id, tranches_id, config):
    """
    Apply variant quality score recalibration for indel variants.

    :param job: Job instance
    :param vcf_id str: vcf file store ID
    :param recal_id: recalibration table file store ID
    :param tranches_id: tranches file store ID
    :param config dict: pipeline configuration options and shared files
    :return str: INDEL recalibrated VCF file store ID
    """
    job.fileStore.logToMaster('Running GATK ApplyRecalibration (INDEL Mode): {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    uuid = config['uuid']
    suffix = config['suffix']

    references = ['genome.fa', 'genome.fa.fai', 'genome.dict']
    inputs = {'toil.vcf': vcf_id,
              'HAPINDEL.recal': recal_id,
              'HAPINDEL.tranches': tranches_id}
    inputs.update({key: config[key] for key in references})
    get_files_from_filestore(job, work_dir, inputs)

    command = ['-T', 'ApplyRecalibration',
               '-input', 'toil.vcf',
               '-o', 'toil.vqsr.vcf',
               '-R', 'genome.fa',
               '-nt', '1',
               '-ts_filter_level', '99.0',
               '-tranchesFile', 'HAPINDEL.tranches',
               '-recalFile', 'HAPINDEL.recal',
               '-mode', 'INDEL']

    if config['unsafe_mode']:
        command = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY'] + command

    outputs={'toil.vqsr.vcf': None}
    docker_call(work_dir = work_dir,
                env={'_JAVA_OPTIONS':'-Djava.io.tmpdir=/data/ -Xmx{}'.format(config['xmx'])},
                parameters = command,
                tool = 'quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs = inputs.keys(),
                outputs = outputs,
                mock=config['mock_mode'])
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'toil.vqsr.vcf'))


