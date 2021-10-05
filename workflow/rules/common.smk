# Python functions used throughout snakemake workflow

def get_raw_reads(wildcards):
    """
    Extract forward and reverse read FASTQ paths from file
    """
    if wildcards.sample.startswith('s_'):
        R1 = glob.glob(TOR_RAW_READ_DIR + '/{0}/{0}_*_1.fq.gz'.format(wildcards.sample))[0]
        R2 = glob.glob(TOR_RAW_READ_DIR + '/{0}/{0}_*_2.fq.gz'.format(wildcards.sample))[0]
    else:
        R1 = glob.glob(GLUE_RAW_READ_DIR + '/{0}/{0}_*_1.fq.gz'.format(wildcards.sample))[0]
        R2 = glob.glob(GLUE_RAW_READ_DIR + '/{0}/{0}_*_2.fq.gz'.format(wildcards.sample))[0]
    return { 'read1' : R1, 'read2' : R2 }

def get_toronto_bam(wildcards):
    """
    Returns Toronto BAM for samples to be included in GLUE
    """
    all_bams = expand(rules.samtools_markdup.output.bam, sample = SAMPLES)
    bam = [bam for bam in all_bams if os.path.basename(bam).startswith(wildcards.sample)]
    return bam

def get_all_bams(wildcards):
    """
    Returns list with paths to 500 GLUE bams and 20 Downsampled Toronto Bams
    """
    glue_bams = expand(rules.samtools_markdup.output.bam, sample = SAMPLES)
    glue_bams_noTor = [bam for bam in glue_bams if not os.path.basename(bam).startswith('s_')]
    tor_bams = expand(rules.downsample_toronto_bam.output, sample = TOR_SAMPLES)
    return glue_bams_noTor + tor_bams

