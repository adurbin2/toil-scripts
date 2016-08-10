import textwrap


def generate_config():
    return textwrap.dedent("""
        # GATK Germline Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline
        # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank.
        ##############################################################################################################
        # Required: Reference Genome URL (hg19.fa)
        genome-fasta:

        # Optional: Reference Genome Index URL (hg19.fa.fai)
        genome-fai:

        # Optional: Reference Genome Sequence Dictionary URL (hg19.dict)
        genome-dict:

        # Optional: 1000G INDELs URL (1000G_phase1.indels.hg19.sites.fixed.vcf)
        phase:

        # Optional: Mills INDELs URL (Mills_and_1000G_gold_standard.indels.hg19.sites.vcf)
        mills:

        # Optional: dbSNP URL (dbsnp_138.hg19.excluding_sites_after_129.vcf)
        dbsnp:

        # Optional: HapMap URL (hapmap_3.3.hg19.sites.vcf)
        hapmap:

        # Optional: Omni URL (1000G_omni2.5.hg19.sites.vcf)
        omni:

        # Optional: Run BWA on fastqs (Boolean)
        run-bwa:

        # Optional. BWA trims adapters (Boolean)
        trim: True

        # Optional: BWA Index amb file
        amb:

        # Optional: BWA Index ann file
        ann:

        # Optional: BWA Index bwt file
        bwt:

        # Optional: BWA Index pac file
        pac:

        # Optional: BWA Index sa file
        sa:

        # Optional: BWA Index alt file
        alt:

        # Optional: Run GATK Preprocessing (Boolean)
        preprocess:

        # Optional: Run GATK VQSR (Boolean)
        run-vqsr:

        # Optional: Joint Calling (Boolean)
        joint:

        # Optional: Run Oncotator (Boolean)
        run-oncotator:

        # Optional: Oncotator Database (URL)
        oncotator-db:

        # Optional: Approximate input file size (i.e. 50G)
        file-size: 50G

        # Memory allocation for Java option Xmx (i.e. 100G)
        xmx: 100G

        # Optional: Suffix to be added to output filename (i.e. .toil)
        suffix:

        # Optional: Path to output directory (PATH/URL)
        output-dir:

        # Optional: Synapse username and password
        synapse-name:
        synapse-pwd:

        # Optional: Path to key file for SSE-C Encryption
        ssec:

        # Optional: Allow seq dict incompatibility (Boolean)
        unsafe-mode:
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