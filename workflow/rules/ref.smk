# Rules to handle reference genome and annotation file

rule clone_degeneracy:
    """
    Clone Degeneracy GitHub repo for getting 4fold and 0fold sites.
    """
    output:
        temp(directory('Degeneracy'))
    log: 'logs/clone_degeneracy/clone_degeneracy.log'
    shell:
        """
        git clone https://github.com/James-S-Santangelo/Degeneracy.git
        """

rule get_fourfold_zerofold:
    """
    Uses get_4fold_sites.sh from Degeneracy to get 4fold and 0fold sites across white clover
    genome from reference sequence (FASTA) and annotation file (GFF).
    """
    input:
        degen = rules.clone_degeneracy.output,
        ref = REFERENCE_GENOME,
        gff = GFF_FILE
    output:
        expand('{0}/4fold_0fold/Trepens_{{site}}.bed'.format(PROGRAM_RESOURCE_DIR), site=['0fold','4fold'])
    log: 'logs/4fold_0fold/get_fourfold_zerofold.log'
    conda: '../envs/ref.yaml'
    params:
        outpath = '{0}/4fold_0fold/'.format(PROGRAM_RESOURCE_DIR)
    resources:
        mem_mb = 4000,
        time = '06:00:00'
    shell:
        """
        OUT_ABS=$( cd {params.outpath}; pwd )
        ( cd {input.degen} &&
        bash get_4fold_sites.sh {input.gff} {input.ref} $OUT_ABS &&
        rm -rf Degeneracy ) 2> {log}
        """

