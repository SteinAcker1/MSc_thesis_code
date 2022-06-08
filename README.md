# Introduction

This repository contains the complete code needed to replicate the analysis performed in Stein Acker's M.Sc. thesis. Please note that running this code requires some setup and expects the user to have access to a high-performance computing cluster (HPC). The complete code for this project can be found at https://github.com/SteinAcker1/MSc_thesis_code.

# Workflow

This project consists of two separate phases. They are:

## HPC phase

### Introduction
This phase handles the raw FASTQ output from each experiment and is expected to take place in an HPC environment. It was run with Python v3.6.6. Each of the three experiments has its own directory. Within this directory, there is a `main.py` script, as well as a `config` directory. Within the `config` directory, there are two JSON files: `lunarc.json` and `software.json`. The `lunarc.json` file contains HPC-specific information to be used in the generated SLURM scripts and may be edited to work with whatever hardware you are using to perform this analysis. The `software.json` file, on the other hand, contains the specific softwares and version numbers to be used with each step in the `main.py` script. _It should not be edited if you wish to achieve similar results to those in the thesis_.

### Running the code
When you look at each `main.py` script, you will see a section titled `# Setting up file paths` at the top. These filepaths are safe to change to match the layout of your filesystem better. Specifically, you may want to change:

- `raw`: This is the filepath to the raw data you intend to analyze. (Due to quirks of the trusTEr software, it must be in a Python list format, like `["path/to/files"]` rather than simply `"path/to/files"`.)

- `gene_gtf`: This is a path to a GTF file containing genomic gene annotations. The GTF file used here was the Gencode comprehensive gene annotation for all regions v38, and a link to an FTP server for download can be found [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz)

- `te_gtf`: This is a path to a GTF file containing genomic TE annotations. A link to the GTF file used can be found [here](https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_GENCODE_rmsk_TE.gtf.gz)

- `star_index`: This is a path to the STAR index you wish to use in your analysis. This index was created in STAR v2.7.8a using the human genome assembly 38 FASTA file (found [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz) and the `gene_gtf` annotation file listed above.

- `cr_index`: This is a path to a Cell Ranger index.

After you have run the trusTEr workflow, you can download the RData files found in the `/3_mergeSamples/` directory. There will be three outputted RData files: one only containing gene expression information, named with the pattern (experiment).rds (for instance, `parkinsons.rds`) and two RData files with TE data appended, named with the pattern (experiment)\_(condition).rds (for instance, `parkinsons_Ctl.rds`). Download all three files for each experiment.

## Local phase

This phase performs statistical analysis on the RData files produced in the HPC phase. All analysis was performed in R v4.1.2.

### Installing required packages

You can install all required packages in the correct versions using the `easy_install.R` script, though you may of course install the packages by whatever means you like. The required packages are:

- Seurat v4.1.0

- tidyverse v1.3.1

- ggpubr v0.4.0

- viridis v0.6.2

- clusterProfiler v4.2.2

- org.Hs.eg.db v3.14.0

Most of the required packages may also be installed in the correct version in Conda with this command:
```
conda create -n acker_thesis -c conda-forge -c bioconda r-base=4.1.2 r-seurat=4.1.0 r-tidyverse=1.3.1 r-ggpubr=0.4.0 r-viridis=0.6.2 bioconductor-org.Hs.eg.db=3.14.0
```

Unfortunately, clusterProfiler v4.2.2 is not yet available in the Bioconda repository at the time of writing, so you may need to install it on your own. Additionally, a version of Seurat compatible with Apple Silicon Macs is not available in the Conda-Forge repository yet either, so you may need to install that package manually as well if using Apple Silicon.

### Using the code

Some of the filepaths used in the code are specific to my machine. You will need to edit them for the code to run properly. The variables which must be changed are the following:

- `outdir`: Must be changed to a directory on your computer where you want the output deposited

- `parkinsons_sn`, `parkinsons_ctx`, `alzheimers`, and `sample[[x]]`: Must be changed to the location where you saved the RData files from the previous phase, surrounded by `readRDS()`. (Note: For `alzheimers_TEs`, the `files` variable plays a similar role and should be defined as a string filepath that leads to a directory containing both `alzheimers_AD.rds` and `alzheimers_NC.rds`.)

All figures will be found in the location you defined as `outdir`.

The original analysis was performed using RMarkdown documents; however, for your convenience, I have also provided equivalent R scripts created using `knitr::purl()` so that these analyses may be run from a command line.
