rule download_reference:
    """
    Download white clover reference genome from NCBI's FTP directory
    """
    output:
        '{0}/GCA_005869975.1_AgR_To_v5_genomic.fna.gz'.format(REF_DIR)
    log: 'logs/download_reference/download_reference.log'
    params:
        outpath = '{0}/'.format(REF_DIR)
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/005/869/975/GCA_005869975.1_AgR_To_v5/GCA_005869975.1_AgR_To_v5_genomic.fna.gz -P {params.outpath} 2> {log}
        """

rule unzip_reference:
    """
    Unzip reference genome
    """
    input:
        rules.download_reference.output
    output:
        '{0}/GCA_005869975.1_AgR_To_v5_genomic.fna'.format(REF_DIR)
    log: 'logs/unzip_reference/unzip_reference.log'
    shell:
        """
        gunzip {input} 2> {log}
        """

rule samtools_index_reference:
    """
    Index reference genome with samtools
    """
    input:
        rules.unzip_reference.output
    output:
        '{0}/GCA_005869975.1_AgR_To_v5_genomic.fna.fai'.format(REF_DIR)
    conda: '../envs/ref.yaml'
    log: 'logs/samtools_index_reference/index_reference.log'
    shell:
        """
        samtools faidx {input} 2> {log}
        """

rule bwa_index_ref:
    """
    Index reference with bwa to get ready for read mapping
    """
    input:
        rules.unzip_reference.output
    output:
        multiext('{0}/GCA_005869975.1_AgR_To_v5_genomic.fna'.format(REF_DIR), '.amb', '.ann', '.bwt', '.pac', '.sa')
    conda: '../envs/ref.yaml'
    log: 'logs/bwa_index_ref/bwa_index_ref.log'
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    shell:
        """
        bwa index {input} 2> {log}
        """

rule ref_done:
    """
    Writes empty flag file to signal successful completion of reference genome rules.
    """
    input:
        rules.samtools_index_reference.output,
        rules.bwa_index_ref.output
    output:
        '{0}/ref.done'.format(REF_DIR)
    shell:
        """
        touch {output}
        """
