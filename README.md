A Snakemake pipeline to perform some basic QC of raw, trimmed, and mapped Illumina reads for the Global Urban Evolution Project (GLUE)
This pipeline is is automatically pulled by other projects using Snakemake's `module` functionality introduced in version Snakemake v6.

### Overview of the pipeline

The QC pipeline performs the following steps:

1. QC of raw reads using `FastQC`
2. Trimming of raw reads using `fastp`
3. Map reads to white clover reference using `bwa`
4. Sort, index, and mark duplicates in BAM files using `samtools`
5. QC of mapped reads using `bamtools`, `Qualimap`, `bamUtil`, and `multiQC`

### Overview of directories

- [config](./config): Snakemake configuration files for different clusters and programs (e.g., `multiQC`)
- [resources](./resources): Text files used in pipeline (e.g., sample information, chromosomes, etc.)
- [workflow](./workflow): Main Snakemake workflow with rules, environments, scripts, notebooks, and cluster profiles for running the pipeline on Compute Canada SLURM-based clusters.

### Using the pipeline

This pipeline requires `Conda`. A minimal installation of `Conda` (i.e., Miniconda) can be installed by following the instructions for your platform [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Assuming `Conda` is installed, this repository's `Conda` environment can be replicated by running the following command:

`conda env create -f environment.yaml -n glue`

This will create a `Conda` environment named _glue_ containing a minimal set of dependencies required to run the pipeline (e.g., Python 3.8.6 and Snakemake 6.9.1).

After activating the environment (`conda activate glue`), the pipeline can be executed from the [workflow](./workflow) directory by running a command that looks something like:

`snakemake --use-conda --configfile ../config/<configfile> --notemp -j <cores>`

for local execution. Here, `<configfile>` is one of the configfiles in the [config](./config) directory that needs to be modified to match the paths on your system, and `<cores>` is the number of cores available for executing parallel processes. 

For execution on a SLURM cluster, the pipeline can be executed by running:

`snakemake --profile compute-canada --configfile ../config/<configfile>`
