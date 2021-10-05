rule fastqc_raw_reads:
    input:
        unpack(get_raw_reads),
        tmp = rules.create_tmp_dir.output
    output:
        html1 = '{0}/fastqc_raw_reads/{{sample}}_1_fastqc.html'.format(QC_DIR),
        html2 = '{0}/fastqc_raw_reads/{{sample}}_2_fastqc.html'.format(QC_DIR),
        zip1 = temp('{0}/fastqc_raw_reads/{{sample}}_1_fastqc.zip'.format(QC_DIR)),
        zip2 = temp('{0}/fastqc_raw_reads/{{sample}}_2_fastqc.zip'.format(QC_DIR))
    conda: '../envs/qc.yaml'
    log: 'logs/fastqc_raw_reads/{sample}_fastqc_raw_reads.log'
    threads: 2
    resources: 
        mem_mb = 1000, 
        time = '01:00:00'
    shell:
        """
        ( fastqc --threads {{threads}} --outdir {0}/fastqc_raw_reads --noextract --quiet --dir {{input.tmp}} {{input.read1}} {{input.read2}} &&
        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_1_fastqc.html {0}/fastqc_raw_reads/{{wildcards.sample}}_1_fastqc.html &&
        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_1_fastqc.zip {0}/fastqc_raw_reads/{{wildcards.sample}}_1_fastqc.zip &&
        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_2_fastqc.html {0}/fastqc_raw_reads/{{wildcards.sample}}_2_fastqc.html &&
        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_2_fastqc.zip {0}/fastqc_raw_reads/{{wildcards.sample}}_2_fastqc.zip ) 2> {{log}}
        """.format(QC_DIR)

rule fastqc_trimmed_reads:
    input:
        tmp = rules.create_tmp_dir.output,
        read1 = rules.fastp_trim.output.r1_trim,
        read2 = rules.fastp_trim.output.r2_trim
    output:
        html1 = '{0}/fastqc_trimmed_reads/{{sample}}_trimmed_1_fastqc.html'.format(QC_DIR),
        html2 = '{0}/fastqc_trimmed_reads/{{sample}}_trimmed_2_fastqc.html'.format(QC_DIR),
        zip1 = temp('{0}/fastqc_trimmed_reads/{{sample}}_trimmed_1_fastqc.zip'.format(QC_DIR)),
        zip2 = temp('{0}/fastqc_trimmed_reads/{{sample}}_trimmed_2_fastqc.zip'.format(QC_DIR))
    conda: '../envs/qc.yaml'
    log: 'logs/fastqc_trimmed_reads/{sample}_fastqc_trimmed_reads.log'
    threads: 2
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        fastqc --threads {{threads}} --outdir {0}/fastqc_trimmed_reads --noextract --quiet --dir {1} {{input.read1}} {{input.read2}} 2> {{log}}
        """.format(QC_DIR, TMPDIR)

rule qualimap_bam_qc:
    input:
        bam = rules.samtools_markdup.output.bam,
        index = rules.index_bam.output
    output:
        directory('{0}/qualimap/{{sample}}_qualimap_bamqc'.format(QC_DIR))
    log: 'logs/qualimap/{sample}_bamqc.log'
    conda: '../envs/qc.yaml'
    threads: 8
    resources:
        mem_mb = lambda wildcards, threads, input, attempt: attempt * 10000,
        time = lambda wildcards, attempt: str(attempt * 2) + ":00:00"
    shell:
        """
        unset DISPLAY;
        qualimap bamqc -bam {{input.bam}} \
            --paint-chromosome-limits \
            --collect-overlap-pairs \
            -nt {{threads}} \
            -outdir {0}/qualimap/{{wildcards.sample}}_qualimap_bamqc \
            -outformat html \
            --java-mem-size={{resources.mem_mb}}M >> {{log}} 2>&1
        """.format(QC_DIR)

rule bamtools_stats:
    input:
        bam = rules.samtools_markdup.output.bam,
        index = rules.index_bam.output
    output:
        '{0}/bamtools_stats/{{sample}}_bamtools.stats'.format(QC_DIR)
    conda: '../envs/qc.yaml'
    log: 'logs/bamtools_stats/{sample}_bamtools_stats.log'
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    shell:
        """
        bamtools stats -in {input.bam} > {output} 2> {log}
        """

rule bamutil_validate:
    input:
        bam = rules.samtools_markdup.output.bam,
        index = rules.index_bam.output
    output:
        '{0}/bamutil_validate/{{sample}}_validation.txt'.format(QC_DIR)
    log: 'logs/bamutil_validate/{sample}_validation.log'
    conda: '../envs/qc.yaml'
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    shell:
        """
        bam validate --in {input.bam} \
            --so_coord \
            --verbose 2> {output}
        """

rule multiqc:
    """
    Generate single HTML report with all QC info for all samples using multiQC.
    Inputs only enforce dependencies. MultiQC takes results folder as input (i.e., QC_DIR)
    """
    input:
       fastqc_raw = expand(rules.fastqc_raw_reads.output.zip1, sample=SAMPLES),
       fastqc_trim = expand(rules.fastqc_raw_reads.output.zip1, sample=SAMPLES),
       fastp = expand(rules.fastp_trim.output.json, sample=SAMPLES),
       qualimap = expand(rules.qualimap_bam_qc.output, sample=SAMPLES),
       bamstats = expand(rules.bamtools_stats.output, sample=SAMPLES),
       bamutil = expand(rules.bamutil_validate.output, sample=SAMPLES)
    output:
        '{0}/multiqc/multiqc_report.html'.format(QC_DIR)
    conda: '../envs/qc.yaml'
    log: 'logs/multiqc/multiqc.log'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        multiqc --verbose \
            --dirs \
            --force \
            --outdir {0}/multiqc \
            --config ../config/multiqc_config.yaml \
            {0} 2> {{log}}
        """.format(QC_DIR)
