# WES GATK pipeline version 2

## Description
This is a Nextflow based workflow for processing whole exome sequencing (WES) data from parent-child trios (and quads) following GATK4 best practices. The first part of this pipeline performs alignment to the reference genome through HaplotypeCaller. The second part of the pipeline performs variant filtering to get de novo SNPs, de novo indels, homozygous recessive, compound heterozygous, X-linked recessive, and dominant variants.

This pipeline supports running locally, in an interactive session on a HPC, and headless on a HPC. Local use is primarily for testing or running a small number of samples. HPC use is recommended for large batches of samples as the pipeline takes advantage of SLURM job scheduling and can run many processes in parallel. A Docker image is provided which ensures the environment is consistent and portable. The Docker image is pulled into Singularity for running on HPCs.

## Pipeline overview
- Index reference genome (`bwa index`) (Optional)
- Map reads to reference genome (`bwa mem`)
- Mark duplicate reads (`picard MarkDuplicates`)
- Base quality score recalibration (BQSR) (`gatk BaseRecalibrator`)
- HaplotypeCaller (`gatk HaplotypeCaller`)
- Combine trios into GVCF format and call genotypes (`gatk CombineGVCFs` and `gatk GenotypeGVCFs`)
- Variant filtering for de novo SNPs, de novo indels, homozygous recessive, compound heterozygous, X-linked recessive, and dominant variants. Uses ANNOVAR for annotation and BRAVO for allele frequencies.

## Quick Start
1. Create conda environment from included `environment.yml` file and activate. This just gives us access to the `nextflow` command. Activating the environment is not required for the `submitNextflowRun.sh` command since the script already performs it.
- For Agave:
    ```bash
    conda env create -f environment.yml
    source activate wes_gatk_v2
    ```
- For Sol:
    ```bash
    mamba env create -f environment.yml
    source activate wes_gatk_v2
    ```
2. Install or activate Docker or Singularity
3. Download this repository
4. Edit `samplesheet.tsv` with correct sample directory names
5. Point `nextflow.config` file to correct `resourcesDir`, `inDataDir`, `samplesheet`, `pedigreeDir`, and `referenceGenome`
    - Resources required are available on Agave under `/data/kruerlab/refs/wes_gatk_v2/`
    - May have to change `intervalList` depending on kit used for sequencing
6. If running locally, edit `nextflow.config` `standard` profile with appropriate `process.cpus` and `process.memory`
7. Run pipeline
    - If running locally: `nextflow run main.nf -profile local`
    - If running in HPC interactive session: `nextflow run main.nf -profile clusterLocal`
        - Ensure Singularity is loaded `module load singularity/3.8.0`
    - If running headless on HPC (recommended): 
        - Agave: `sbatch submitNextflowRun.sh`
        - Sol: `sbatch submitNextflowRunSol.sh`

## Sample sheet format
The sample sheet is a .tsv file containing a trio on each row. The sample sheet also supports larger families with the other members specified under the `Other_Members` column separated by commas. Each entry contains the directory name that the read pairs are inside of.

Example directory structure:
```
└── data/
    ├── M_F279-001-U/
    │   ├── M_F279-001-U_AH3LJLDSXX_L004_R1_001.fastq.gz
    │   └── M_F279-001-U_AH3LJLDSXX_L004_R2_001.fastq.gz
    ├── M_F279-002-U/
    │   ├── M_F279-002-U_AH3LJLDSXX_L004_R1_001.fastq.gz
    │   └── M_F279-002-U_AH3LJLDSXX_L004_R2_001.fastq.gz
    └── M_F279-003-A/
        ├── M_F279-003-A_AH3LJLDSXX_L004_R1_001.fastq.gz
        └── M_F279-003-A_AH3LJLDSXX_L004_R2_001.fastq.gz
```
The `inDataDir` in the `nextflow.config` file would point to the `data/` directory. `samplesheet.tsv` would contain:

|Mother         | Father      | Child_Affected | Other_Members |
| --------------|-------------|----------------|---------------|
| M_F279-001-U | M_F279-002-U | M_F279-003-A   | M_F279-004-U,M_F279-005-A

## Pipeline outputs
Output directories are specified in the `nextflow.config` file under `outDataDir` and `outTrioDir`. The `outDataDir` contains primarily BAM and SAM files while the `outTrioDir` contains the variant call sheets.

## Other notes
Since this pipeline runs off Nextflow, it comes with some nice functionalities such as the ability to resume runs. This is done by adding the `-resume` flag after the `nextflow run` command. This is useful in cases where the pipeline ran out of time on the cluster or if you forgot some sort of input file (though I did my best on input file checking before processing). If you need to do this, make sure to not delete the `work` directory that is created. It is a temporary directory created by the pipeline and acts as a cache to allow for the resuming of the workflow.

The `nextflow.config` file is very flexible and it is easy to adjust directory paths and resource allocations. Do not feel limited by the default directory locations. These merely act as a starting point for working with this workflow.

Slurm jobs process settings are specified in the `nextflow.config` file under the `cluster` scope. These may be altered if needed if jobs are failing due to insufficient resource allocations. Jobs are automatically resubmitted with more time allowances which should cover edge cases such as large samples taking an unusually large amount of time.

## Manual resource installation
1. Install the `gsutil` tool: https://cloud.google.com/storage/docs/gsutil_install
2. Download reference genome (hg38)
    - `gsutil -m cp -r gs://gcp-public-data--broad-references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy/v0 .`
    - Also available at: https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy
3. Index the reference genome
    - `bwa index Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta`
4. Download dbSNP
    - `gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf`
    - `gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx`
    - Also available at: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0
5. Download ANNOVAR datasets (Warning >700GB)
    ```
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
    perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar esp6500siv2_all  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad211_exome  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20220320  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ljb26_all  humandb/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20180118  humandb/
    ```
    
## Future additions
- Run pipeline on all samples in specified directory
- Automatic concatenation of resequenced samples
