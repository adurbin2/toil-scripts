import textwrap

def generate_config():
    return textwrap.dedent("""
        # GATK Germline Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Reference Genome URL (hg19.fa)
        genome.fa:

        # Optional: Reference Genome Index URL (hg19.fa.fai)
        genome.fa.fai:

        # Optional: Reference Genome Sequence Dictionary (hg19.dict)
        genome.dict:

        # Required: URL (1000G_phase1.indels.hg19.sites.fixed.vcf)
        phase.vcf:

        # Required: URL (Mills_and_1000G_gold_standard.indels.hg19.sites.vcf)
        mills.vcf:

        # Required: URL (dbsnp_132_b37.leftAligned.vcf URL)
        dbsnp.vcf:

        # Required: URL (hapmap_3.3.b37.vcf)
        hapmap.vcf:

        # Required: URL (1000G_omni.5.b37.vcf)
        omni.vcf:

        # Optional: (boolean) Run BWA on fastqs
        run_bwa:

        # Optional: (boolean) Run GATK Preprocessing
        preprocess:

        # Optional: (boolean) Run GATK VQSR
        run_vqsr:

        # Optional: If true, sorts bam
        sort: True

        # Optional. If true, trims adapters
        trim: false

        # Optional: Reference fasta file (amb) -- if not present will be generated
        amb: s3://cgl-pipeline-inputs/alignment/hg19.fa.amb

        # Optional: Reference fasta file (ann) -- If not present will be generated
        ann: s3://cgl-pipeline-inputs/alignment/hg19.fa.ann

        # Optional: Reference fasta file (bwt) -- If not present will be generated
        bwt: s3://cgl-pipeline-inputs/alignment/hg19.fa.bwt

        # Optional: Reference fasta file (pac) -- If not present will be generated
        pac: s3://cgl-pipeline-inputs/alignment/hg19.fa.pac

        # Optional: Reference fasta file (sa) -- If not present will be generated
        sa: s3://cgl-pipeline-inputs/alignment/hg19.fa.sa

        # Optional: Reference fasta file (fai) -- If not present will be generated
        fai: s3://cgl-pipeline-inputs/alignment/hg19.fa.fai

        # Optional: (string) Path to Key File for SSE-C Encryption
        ssec:

        # Optional: Alternate file for reference build (alt). Necessary for alt aware alignment
        alt:

        # Optional: Synapse account and password
        synapse_id:
        synapse_pwd:


        # Required: Approximate input file size. Provided as a number followed by (base-10) [TGMK]. E.g. 10M, 150G
        file-size: 50G

        # Memory allocation for Java option Xmx
        xmx: 100G

        # Optional: (string) Suffix to be added to final output
        suffix:

        # Optional: (string) Path to output directory
        output_dir:

        # Optional: (boolean) Set to True to allow seq dict incompatibility
        unsafe_mode:

        # Optional: (boolean) Set to True to run pipeline in mock mode
        mock_mode:
        """[1:])

def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #   There are 2 tab-separated columns: UUID and URL
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   URL         A URL (http://, ftp://, file://, s3://) pointing to the input SAM or BAM file
        #
        #   Example below. Lines beginning with # are ignored.
        #
        #   UUID_1    file:///path/to/sample.bam
        #                   OR
        #   UUID_1    file:///path/to/sample.1.fq @RG\tID:foo\tSM:bar
        #
        #   Place your samples below, one per line.
    """[1:])