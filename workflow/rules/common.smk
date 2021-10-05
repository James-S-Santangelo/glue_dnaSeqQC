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

def get_fastas_to_concat(wildcards):
    """
    Retrieves all consensus FASTA files for rbcL or matK, depending on value of "gene" wildcard.
    """
    if wildcards.gene == 'rbcl':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='rbcl')
    elif wildcards.gene == 'matk':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='matk')

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

def get_bed_to_subset(wildcards):
    """
    Returns BED file for either 0fold or 4fold sites, depending on value of "site" wildcard.
    """
    all_bed_files = rules.get_fourfold_zerofold.output
    bed = [bed for bed in all_bed_files if wildcards.site in os.path.basename(bed)]
    return bed

def get_bams_for_angsd(wildcards):
    """
    Returns correct text file with paths to BAM file to be including in ANGSD, depending on
    value of "sample_set" wildcard.
    """
    if wildcards.sample_set == 'highErrorRemoved':
        return rules.create_bam_list_highErrorRemoved.output
    elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
        return rules.create_bam_list_finalSamples_lowCovRemoved.output

def angsd_sfs_input(wildcards):
    """
    Returns dictionary with correct SAF and sites files, depending 
    on whether all sites are being analysed (i.e., "site" wildcard == "allSites")
    or whether degenerate sites are being analysed ("site" wildcard == '4fold' or '0fold')
    """
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output.sites
        idx = rules.angsd_index_allSites.output.idx
    else:
        sites = rules.split_angsd_sites_byChrom.output.sites
        idx = rules.angsd_index_degenerate.output.idx
    return { 'saf_idx' : saf_idx, 'sites' : sites, 'sites_idx' : idx }

def angsd_estimate_thetas_input(wildcards):
    """
    Returns dictionary with correct SAF, sites, and sfs files, depending 
    on whether all sites are being analysed (i.e., "site" wildcard == "allSites")
    or whether degenerate sites are being analysed ("site" wildcard == '4fold' or '0fold')
    """
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    sfs = rules.angsd_estimate_sfs.output
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output.sites
        idx = rules.angsd_index_allSites.output.idx
    else:
        sites = rules.split_angsd_sites_byChrom.output.sites
        idx = rules.angsd_index_degenerate.output.idx
    return { 'saf_idx' : saf_idx, 'sfs' : sfs, 'sites' : sites, 'sites_idx' : idx }

def get_angsd_stats_toConcat(wildcards):
    """
    Returns list with correct diversity and neutrality stats files for concatenation, depending on combination
    of "sample_set" and "site" wildcard values
    """
    if wildcards.sample_set == 'highErrorRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='0fold', sample_set='highErrorRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='4fold', sample_set='highErrorRemoved')
        else:
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='allSites', sample_set='highErrorRemoved')
    elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_lowCovRemoved')
        else:
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_lowCovRemoved')

def get_angsd_sfs_toConcat(wildcards):
    """
    Returns list with correct sfs files for concatenation, depending on combination
    of "sample_set" and "site" wildcard values
    """
    if wildcards.sample_set == 'highErrorRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='0fold', sample_set='highErrorRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='4fold', sample_set='highErrorRemoved')
        else:
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='allSites', sample_set='highErrorRemoved')
    elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_lowCovRemoved')
        else:
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_lowCovRemoved')

def get_angsd_gl_toConcat(wildcards):
    """
    Returns list with correct genotype likelihood files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    if wildcards.maf == '0.005':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.01':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.05':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    return out

def get_angsd_maf_toConcat(wildcards):
    """
    Returns list with correct minor allele frequency files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    if wildcards.maf == '0.005':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.01':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.05':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    return out

def get_habitat_saf_files_byCity(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx, city=CITIES, habitat=HABITATS, site=['4fold'])
    city_saf_files = [x for x in all_saf_files if wildcards.city in x and wildcards.site in x]
    return city_saf_files

def get_bamLists_toConcat(wildcards):
    """
    Collect text files with paths to urban and rural bams by city
    """
    all_bam_lists = expand(rules.create_bam_list_byCity_byHabitat.output, city = wildcards.city, habitat = ['u', 'r'])
    return all_bam_lists

def get_files_saf_estimation_byPopulation(wildcards):
    """
    Collect files required to estimate SAF likelihoods by population in each city. Population ID extracted
    from files created by Checkpoint.
    """
    bams = '{0}/bam_lists/by_city/{{city}}/by_pop/{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019101.1', site='4fold')
    sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019101.1', site='4fold') 
    ref = REFERENCE_GENOME
    return { 'bams' : bams, 'sites_idx' : sites_idx, 'sites' : sites, 'ref' : ref }
    
def aggregate_input_theta(wildcards):
    """
    Collect population ID ('popu') wildcard values from checkpoint.
    Collects output for estimation of thetas
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    return expand('{0}/summary_stats/thetas/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR), city=wildcards.city, popu=pops, site='4fold')

def get_population_saf_files_byCity(wildcards):
    """
    Get SAF files for two populations for which to estimate joint SFS.
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byPopulation.output.saf_idx, city = wildcards.city, popu=pops, site='4fold')
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '_{0}_'.format(pop1) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '_{0}_'.format(pop2) in os.path.basename(x)]
    return saf1 + saf2

def get_urban_rural_bam_lists(wildcards):
    """
    Collect files with paths to urban and rural bams by City. Return as dictionary. 
    """
    urban = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat='u')[0]
    rural = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat='r')[0]
    return { 'urban_bams' : urban, 'rural_bams' : rural }

def get_population_saf_and_sfs_files_byCity(wildcards):
    """
    Get SAF and SFS files for two populations for which Fst should be estimated. 
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byPopulation.output.saf_idx, city = wildcards.city, popu=pops, site='4fold')
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '_{0}_'.format(pop1) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '_{0}_'.format(pop2) in os.path.basename(x)]
    saf_files = saf1 + saf2
    sfs = expand(rules.angsd_estimate_joint_sfs_populations.output, city = wildcards.city, site='4fold', pop_comb=wildcards.pop_comb)
    return { 'saf_files' : saf_files, 'sfs' : sfs }

def get_files_for_saf_estimation_byHabitat(wildcards):
    """
    Get files to estimate SAF likelihhods for urban and rural habitats by city.
    """
    sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019101.1', site='4fold')
    sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019101.1', site='4fold')
    ref = REFERENCE_GENOME
    bams = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat=wildcards.habitat)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_files_for_permuted_saf_estimation(wildcards):
    """
    Get files to estimate SAF likelihoods for permuted versions of "urban" and "rural" populations
    """
    sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019101.1', site='4fold')
    sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019101.1', site='4fold')
    ref = REFERENCE_GENOME
    if wildcards.habitat == 'u':
        bams = expand(rules.create_random_bam_list_byCity_byHabitat.output.urban, city=wildcards.city, seed=wildcards.seed)
    elif wildcards.habitat == 'r':
        bams = expand(rules.create_random_bam_list_byCity_byHabitat.output.rural, city=wildcards.city, seed=wildcards.seed)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity_permuted(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    all_saf_files = expand(rules.angsd_permuted_saf_likelihood_byCity_byHabitat.output.saf_idx, city=CITIES, habitat=HABITATS, site=['4fold'], seed=wildcards.seed)
    city_saf_files = [x for x in all_saf_files if wildcards.city in x and wildcards.site in x]
    return city_saf_files

def aggregate_input_fst(wildcards):
    """
    Collect population ID ('popu') wildcard values from checkpoint.
    Collects output for estimation of fst
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    pop_combinations = [c[0] + '_' + c[1] for c in list(itertools.combinations(pops, 2))]
    return expand('{0}/summary_stats/fst/fst1/{{city}}/pairwise/{{city}}_{{site}}_{{pop_comb}}_readable.fst'.format(ANGSD_DIR), city=wildcards.city, pop_comb=pop_combinations, site='4fold')

def get_ngsadmix_logfiles_byCity(wildcards):
    """
    Get NGSadmix logfiles for each city. Used to generate input file for CLUMPAK
    """
    return expand(rules.ngsadmix.output.lf, city=wildcards.city, k=NGSADMIX_K, seed=NGSADMIX_SEEDS, site='4fold', maf='0.05')

def get_files_for_saf_estimation_snps_hcn_chroms(wildcards):
    """
    Get files to estimate SAF likelihhods for urban and rural habitats by city.
    """
    if wildcards.gene == 'li':
        sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019108.1', site='4fold')
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019108.1', site='4fold')
    elif wildcards.gene == 'ac':
        sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019103.1', site='4fold')
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019103.1', site='4fold')
    ref = REFERENCE_GENOME
    bams = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat=wildcards.habitat)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity_hcn_chroms(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    saf_files = expand(rules.angsd_saf_likelihood_snps_hcn_chroms.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=['4fold'], gene=wildcards.gene)
    return saf_files
