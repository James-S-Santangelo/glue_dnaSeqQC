# Python functions used throughout snakemake workflow

def run_fast_scandir(dir, sample):
    subfolders, files = [], []

    for f in os.scandir(dir):
        if f.is_dir():
            subfolders.append(f.path)
        if f.is_file():
            if re.findall(r'^{0}_'.format(sample), f.name):
                files.append(f.path)
    for dir in list(subfolders):
        sf, f = run_fast_scandir(dir, sample)
        subfolders.extend(sf)
        files.extend(f) 
    return subfolders, sorted(files)

def get_raw_reads(wildcards):
    folders, reads = run_fast_scandir(RAW_READ_DIR, wildcards.sample)
    R1 = reads[0]
    R2 = reads[1]
    return { 'read1' : R1, 'read2' : R2  }

def get_correct_bam(wildcards):
    """
    Retieve GLUE bam or downsampled Toronto bam based on sample name. Used to pass BAMs to QC
    """
    if wildcards.sample.startswith('s_'):
        bam = expand(rules.downsample_toronto_bam.output, sample=wildcards.sample)
        idx = expand(rules.index_toronto_bam.output, sample=wildcards.sample)
    else:
        bam = expand(rules.samtools_markdup.output.bam, sample=wildcards.sample)
        idx = expand(rules.index_bam.output, sample=wildcards.sample)
    return { 'bam' : bam, 'idx' : idx }

