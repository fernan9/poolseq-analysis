# PoolSeq Variant Calling Pipeline

## Project Overview
This pipeline processes pooled sequencing (PoolSeq) data from raw FastQ files to variant calls (VCF) using SLURM job scheduling. 
The workflow includes quality control, trimming, alignment, sorting, pileup generation, and variant calling.

## High Speed Computation Resources
All data and analysis will be conducted in the LCC super computer of the University of Kentucky via VPN.

### Project directory in nmte222_uksr
This directory has 25TB of space and will hold most data.
`/project/nmte222_uksr/perezgalvez/01.RawData`

However, the project should have the following structure:
<pre><code>/project/nmte222_uksr/perezgalvez/poolseq_project/
├── 1_raw_data/ # Original FastQ files
├── 2_fastqc/ # Quality reports
├── 3_trimmed/ # Trimmed FastQ files
├── 4_aligned/ # Alignment files (sai, sam, bam)
├── 5_sorted/ # Sorted BAM files
├── 6_pileups/ # Pileup files
├── 7_variants/ # Final variant calls
├── 00_scripts/ # SLURM scripts
└── 000_logs/ # SLURM output logs
</code></pre>

## Data transfer
Files were uploaded from Windows SCP. The original directory was organized by sequencing lane, each having an independent file per sequencing sense, was reorganized by extracting all *.fq.gz files into the project directory.

<pre><code>../01_RawData/
├── A0
    ├── A0_H7MY2DSX5_L1_1.fq.gz   # sequencing from one sense
    ├── A0_H7MY2DSX5_L1_2.fq.gz   # sequencing from other sense
</code></pre>
The procedure was made via SHH using the following command that extracts all the names of the Gzip files and uses them as arguments for the move command into the required directory.

`$ find 01.RawData -type f -name "*.gz" -print0 | xargs -0 mv -t poolseq_project/1_raw_data/`


# Pipeline Steps

Programs are called from the Singularity container in LCC:
https://ukyrcd.atlassian.net/wiki/spaces/RCDDocs/pages/162104319/Conda+and+Containers+on+LCC


## Header example

The information of each script is contained in the header, this is an example from the first script which already runs.

<pre><code>#!/bin/bash
#SBATCH --time=08:15:00  
#SBATCH --job-name=01_fastqc
#SBATCH --ntasks=1 
#SBATCH --partition=SKY32M192_L 
#SBATCH --error=000_logs/01_fastqc_%j.err 
#SBATCH --output=000_logs/01_fastqc_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com   </code></pre>

### Partition adjustments
The resources avaiable can be accessed with `$ sinfo`. Select one idle node and use it in `--partition=SKY32M192_L` to adjust the partition to be used in the script.
### Run the script
In order to run a script, use the following format from the working directory:

`$ sbatch first_job.sh`


## 1. Quality Control (FastQC)
This initial step will provide information on the quality of the readings. In general, the sequencing was good in all cases.
- **Script**: `scripts/01_fastqc.sh`
- **Input**: Raw FastQ files in `1_raw_data/`
- **Output**: Quality reports in `2_fastqc/`
- **Tools**: FastQC
- **Parameters**: Default settings


### Pending example of results
Show here an example

## 2. Trimming (PoPoolation v1.2.2 script)

The resequencing readings from the fastQC files must be preprocessed before aligning. This is done with the Popoolation pipeline script `trim-fastq-pl`

- **Script**: `scripts/02_trimming.sh`
- **Input**: Raw FastQ files
- **Output**: Trimmed FastQ files in `3_trimmed/`
- **Tools**: Popoolation trim pearl script
- **Parameters**:
  - --quality-threshold 20
  - --min-length 50
  - --fastq-type sanger

### Script logic

The idea is to process every file in `../1_rawdata` in the same script

## 3. Alignment (BWA)
- **Script**: `scripts/03_alignment.sh`
- **Input**: Trimmed FastQ files
- **Output**: SAM files in `aligned/`
- **Tools**: BWA (aln + sampe)
- **Parameters**:
  - Reference genome: `/path/to/reference_genome.fa`
  - Threads: 8

### 4. SAM to BAM Conversion
- **Script**: `scripts/04_sam_to_bam.sh`
- **Input**: SAM files
- **Output**: BAM files in `aligned/`
- **Tools**: samtools view
- **Parameters**: Default

### 5. Sorting and Indexing
- **Script**: `scripts/05_sort_index.sh`
- **Input**: BAM files
- **Output**: Sorted BAM files and indices in `sorted/`
- **Tools**: samtools sort + index
- **Parameters**: Default

### 6. Pileup Generation
- **Script**: `scripts/06_pileup.sh`
- **Input**: Sorted BAM files
- **Output**: Pileup files in `pileups/`
- **Tools**: samtools mpileup
- **Parameters**:
  - Reference genome required
  - Default quality filters

### 7. Variant Calling (VarScan)
- **Script**: `scripts/07_variant_calling.sh`
- **Input**: Pileup files
- **Output**: VCF files in `variants/`
- **Tools**: VarScan
- **Parameters**:
  - min-coverage: 10
  - min-reads2: 4
  - min-avg-qual: 20
  - p-value: 0.05


text

## Execution Instructions

### Sequential Submission
```bash
# Submit jobs with dependencies
JOB1=$(sbatch scripts/01_fastqc.sh | awk '{print $4}')
JOB2=$(sbatch --dependency=afterok:$JOB1 scripts/02_trimming.sh | awk '{print $4}')
JOB3=$(sbatch --dependency=afterok:$JOB2 scripts/03_alignment.sh | awk '{print $4}')
JOB4=$(sbatch --dependency=afterok:$JOB3 scripts/04_sam_to_bam.sh | awk '{print $4}')
JOB5=$(sbatch --dependency=afterok:$JOB4 scripts/05_sort_index.sh | awk '{print $4}')
JOB6=$(sbatch --dependency=afterok:$JOB5 scripts/06_pileup.sh | awk '{print $4}')
JOB7=$(sbatch --dependency=afterok:$JOB6 scripts/07_variant_calling.sh | awk '{print $4}')
Monitoring Jobs
bash
squeue -u $USER          # View your running jobs
sacct -j <JOBID>         # Check status of specific job
scancel <JOBID>          # Cancel a job
PoolSeq-Specific Considerations
Coverage Requirements: Ensure sufficient coverage across all pools

Allele Frequency: Adjust variant calling thresholds for pooled data

Replicates: Process technical/biological replicates consistently

Population Parameters: Consider pool size when interpreting results

Dependencies
FastQC

Trimmomatic

BWA

SAMtools

VarScan

Version History
2023-11-20: Initial pipeline setup

Notes
Adjust memory/time requirements in SLURM scripts based on your dataset size

Check intermediate files at each step before proceeding

Consider adding steps for coverage analysis if needed

text

This README provides:
1. A clear roadmap of your analysis
2. Documentation of parameters and tools
3. Execution instructions
4. PoolSeq-specific considerations
5. Version tracking

You can extend this with:
- Specific sample information
- Custom parameters you used
- Any troubleshooting notes
- References to important papers or methods