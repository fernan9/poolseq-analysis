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
├── 0_genome/ # Genome file for Index
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

The objective is to call SNP variants from aligned data. Alignment between the trimmed samples and the genome must be computed by pairs and requires previous indexing for the genome. The parameter used to align the resequencing data must be benchmarked to select a balance between precision and speed. Because in pool-seq approahces all the variants come from within the sample, the selection of variants has to be conservative.

- **Script**: `scripts/03_index.sh`
              `scripts/03_alignment.sh`
              `scripts/03_benchmarking.sh`
              `scripts/03_benchmarking_analysis.sh`
- **Input**: Trimmed FastQ files
- **Output**: SAM, BAM, VCF files in `4_aligned/`
- **Tools**: BWA (aln + sampe)
- **Parameters**:
  - Reference genome: `/path/to/reference_genome.fa`
  - Threads: 8

### Reference Genome (Index)
Firts, make an index from the genome reference with bwa index, only once. I used the release 6 from *Drosophila melanogaster* downloaded from flybase.

```bash
# Directories
CONTAINER="/share/singularity/images/ccs/conda/lcc-conda-2-centos8.sinf"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

singularity run --app bwa0717 $CONTAINER bwa index $GENOME
```
### Alignment Parameters

Second, align the trimmed files using bwa aln and then producing the lignment SAM file with sampe, for each  trimmed paired sample. Below is a description of the parameters used in the legacy code.

|Parameters	|Description			     			            |
|-----------|---------------------------------------|
|bwa aln	  |specifies that the "aln" algorithm      |
|-t 8       | number of threads (i.e. eight CPU cores) for the alignment process, potentially speeding up the analysis.|
|-o 2     	|maximum number of gap opens allowed in the alignment process.|
|-d 12     	|maximum gap extension penalty that can be assigned in the alignment process.|
|-e 12     	|maximum edit distance allowed between the read and the reference genome, in terms of the number of mismatches and gaps.|
|-l 100     	|seed length used for the alignment process. In this case, a seed length of 100 base pairs is used.|
|-n 0.01     	|maximum allowed mismatch rate between the read and the reference genome. In this case, 0.01 (or 1%) .|

### Bechmarking
However, we want to know the influence of these parameters in variant calling and potentially biasing the analysis of poolseq. For example, by introducing variation by allowing too many gaps or variants in alignment. Favor specificity (low false positives) over sensitivity.

Th goal with this benchmarking, then, is to minimize false positives in variant calls while detecting low-frequency alleles in pooled samples consistently across replicates.

The alignment was conducted with three parameter sets (permissive, adjusted, strict) to evaluate characteristics on its performance.

|Parameters	|Permissive	  |Adjusted	  |Strict |
|-----------|-------------|-----------|-------|
|-n         |  0.04       |	0.02      |	0.01  |
|-o         |	2	          |1	        |1      |
|-e	        |7	          |5	        |3      |
|-l	        |100	        |32	        |32     |

Comparison between alignment was based on the performance metrics of `samtools flagstat` for alignment statistics, the mean depth distribution and the SNP recovery. The idea is to expect the adjested parameters to recover ≥95% known SNPs in comparison to knwon variants as those in FlyBase (but still a good metric has to be determined, specially because the two experiments are different in their expected genetic diversity).

### Benchmarking results

The sample run was B2a to create the SAM file.


|Parameters	      |Permissive	  |Adjusted	  |Strict |
|-----------------|-------------|-----------|-------|
|Mapped rate      | 93.96%    |	96.22%    |	 96.28%  |
|Singletons       |	1.15%      |	0.92%      | 0.91%   |
|Properly Paired	|	91.74%     |	94.10%   |  94.16%   |
|Read Depth	Mean  | 45.1183    |	46.2965   |   46.324   |
|Read Depth (Std Dev)|(106.188)|	(109.861) |  (110.21) |
|SNPs called	    |	1,197,475  |	1,244,337 | 1,244,103 |

One interesting observation is that the strict parameters calls a larger amount of SNPs as well as an increased mapping rate. One would expect that a permissive set of parameters would allow more SNPs to be read, but instead it appears that the contrary is true, perhaps by increasing ambiguity. This idea is reinforced by the increased amount of properly paired reads in strict and adjusted parameters while as well as the decreased number of singletons in the same parameter sets.

The adjusted and strict parameter sets have similar results, which suggests that the seed length `-l` has a large impact in the alignment results. Perhaps this is also true for the number of gap opens `-o` but the independent effects would be hard to tease appart as a complete mapping of the combinations will be required and at the end the choice of parameters will still be arbitrary, alas based on comparison between results.

For the moment, it seems appropiate te apply the **Adjusted** set.

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