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

<pre><code>../poolseq_project/
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
Note: project was moved to a new location

<pre><code>mv /project/nmte222_uksr/perezgalvez/poolseq_project /mnt/gpfs2_4m/scratch/frpe222/</code></pre>
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

To check the status of the current jobs, consult:

`$ squeue -u frpe222`

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

The fastQ files must be preprocessed before aligning. This is done with the Popoolation pipeline script `trim-fastq.pl`. The equivalent script from the legacy 2023 script is `submit-frpe05.sh` and is located at the **2023_slurm_scripts** directory.

- **Script**: `00_scripts/02_trimming.sh` and  `00_scripts/02_ungzipping.sh`
- **Input**: Raw FastQ files
- **Output**: Trimmed FastQ files in `3_trimmed/`
- **Tools**: Popoolation trim pearl script
- **Parameters**:
  - --quality-threshold 20
  - --min-length 50
  - --fastq-type sanger
  - --output1
  - --output2
  - --outputse

### UnGzipping files
First the *.fq.gz files must be decompressed. The code below was applied as a slurm script named `02_ungzipping.sh`
```
cd /project/nmte222_uksr/perezgalvez/1_raw_data
for f in *.fq.gz; do
  gunzip -c "$f" > "${f%.gz}"
done
```

## 3. Alignment (BWA)

Now, the alignment between the trimmed samples and the genome must be computed by pairs. The legacy script corresponds to XXX

- **Script**: `scripts/03_alignment.sh`
- **Input**: Trimmed FastQ files
- **Output**: SAM files in `aligned/`
- **Tools**: BWA (aln + sampe)
- **Parameters**:
  - Reference genome: `/path/to/reference_genome.fa`
  - Threads: 8

first make an index from the genome reference with bwa index, only once

```console
singularity run --app bwa0717 /share/singularity/images/ccs/conda/amd-conda1-centos8.sinf bwa index /scratch/frpe222/Genome/dmel-short-header.fa
```
second align the trimmed files using bwa sampe, for each  trimmed pairedsample

```python
# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A0_trimed_1 > /scratch/frpe222/Trim/A0_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A0_trimed_2 > /scratch/frpe222/Trim/A0_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A0_trimed_1.sai /scratch/frpe222/Trim/A0_trimed_2.sai /scratch/frpe222/Trim/A0_trimed_1 /scratch/frpe222/Trim/A0_trimed_2 > /scratch/frpe222/Mapped/A0.sam

# "bwa aln" specifies that the "aln" algorithm within the BWA software package should be used for read alignment.
# "-t 8" specifies that the tool should use eight threads (i.e. eight CPU cores) for the alignment process, potentially speeding up the analysis if multiple cores are available.
# "-o 2" specifies the maximum number of gap opens allowed in the alignment process. In this case, up to two gap opens are permitted.
# "-d 12" specifies the maximum gap extension penalty that can be assigned in the alignment process.
# "-e 12" specifies the maximum edit distance allowed between the read and the reference genome, in terms of the number of mismatches and gaps.
# "-l 100" specifies the seed length used for the alignment process. In this case, a seed length of 100 base pairs is used.
# "-n 0.01" specifies the maximum allowed mismatch rate between the read and the reference genome. In this case, a maximum mismatch rate of 0.01 (or 1%) is allowed.
```

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

```bash
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
2025-07-23: Initial pipeline setup

Notes
Adjust memory/time requirements in SLURM scripts based on your dataset size

Check intermediate files at each step before proceeding

Consider adding steps for coverage analysis if needed