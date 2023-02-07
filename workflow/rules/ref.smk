rule samtools_index_reference:
    """
    Index reference genome with samtools
    """
    input:
        REFERENCE_GENOME
    output:
        f'{REFERENCE_GENOME}.fai'
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
        REFERENCE_GENOME
    output:
        multiext(f'{REFERENCE_GENOME}', '.amb', '.ann', '.bwt', '.pac', '.sa')
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
