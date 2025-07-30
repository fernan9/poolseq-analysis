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
- **Output**: SAM files in `4_aligned/`
- **Tools**: BWA (aln + sampe)
- **Parameters**:
  - Reference genome: `/path/to/reference_genome.fa`
  - Threads: 8

first make an index from the genome reference with bwa index, only once

```console
singularity run --app bwa0717 /share/singularity/images/ccs/conda/amd-conda1-centos8.sinf bwa index /scratch/frpe222/Genome/dmel-short-header.fa
```
second align the trimmed files using bwa sampe, for each  trimmed pairedsample. In the case of Poolseq data, cetrtain parameters must be benchmarked to optimize variant call. This will be made by doing a benchmark for two file sets, one from B2 and one from P3.

PoolSeq Alignment Benchmarking Guide

Optimizing bwa aln Parameters for Illumina Resequencing Data

1. Goals
Accuracy: Minimize false positives in variant calls.

Sensitivity: Detect low-frequency alleles in pooled samples.

Reproducibility: Consistent performance across replicates.

2. Test Dataset Preparation
```bash
# Extract 1M reads (250K pairs) for quick testing  
zcat ${INPUT_DIR}/sample_R1.trimmed.fq.gz | head -n 1000000 > test_R1.fq  
zcat ${INPUT_DIR}/sample_R2.trimmed.fq.gz | head -n 1000000 > test_R2.fq  
```
3. Parameter Comparison

A. Run Alignments
```bash
# Strict (conservative)  
singularity exec $CONTAINER bwa aln -n 0.01 -o 1 -e 3 -l 32 $GENOME test_R1.fq > strict_R1.sai  

# Adjusted (recommended for PoolSeq)  
singularity exec $CONTAINER bwa aln -n 0.02 -o 1 -e 5 -l 32 $GENOME test_R1.fq > adjusted_R1.sai  

# Permissive (default-like)  
singularity exec $CONTAINER bwa aln -n 0.04 -o 2 -e 7 -l 100 $GENOME test_R1.fq > permissive_R1.sai  
```
B. Generate SAM Files
```bash
for prefix in strict adjusted permissive; do  
  singularity exec $CONTAINER bwa sampe $GENOME ${prefix}_R1.sai ${prefix}_R2.sai \  
    test_R1.fq test_R2.fq > ${prefix}.sam  
done  
```
4. Key Metrics

A. Alignment Statistics
```bash
for sam in *.sam; do  
  echo "===== $sam ====="  
  samtools flagstat $sam  
  echo "Mapped reads: $(grep -c '^[^@]' $sam)"  
done  
Target: >90% mapping rate.
```
B. Depth Distribution
```bash
for sam in *.sam; do  
  samtools view -b $sam | samtools depth - | awk '{sum+=$3} END {print "Mean depth:", sum/NR}'  
done  
```
C. SNP Recovery (vs. Known Variants)
```bash
for sam in *.sam; do  
  bcftools mpileup -f $GENOME $sam | bcftools call -mv -Ov -o ${sam%.sam}.vcf  
  echo "SNPs in ${sam%.sam}.vcf: $(grep -vc '^#' ${sam%.sam}.vcf)"  
done  
```
Expected: Adjusted parameters recover ≥95% known SNPs (e.g., FlyBase variants).

5. Visual Comparison


Example Output Table

|Parameters	|Mapped Rate	|Mean Depth	|SNPs Detected|
|-----------|-------------|-----------|-------------|
|Strict	    |88%          |	45X       |	9,201       |
|Adjusted   |	93%	        |48X	      |9,812        |
|Permissive	|95%	        |49X	      |10,305 (12% FPs)|


Plotting (R/Python)
```r
# R example (ggplot2)  
library(ggplot2)  
data <- read.table("stats.tsv", header=TRUE)  
ggplot(data, aes(x=Parameters, y=MappedRate)) + geom_bar(stat="identity")
```  
6. Automation Script
Save as benchmark.sh:

```bash
#!/bin/bash  
#SBATCH --job-name=alignment_benchmark  
#SBATCH --time=2:00:00  
#SBATCH --cpus-per-task=8  

declare -A PARAMS=(  
  ["strict"]="-n 0.01 -o 1 -e 3 -l 32"  
  ["adjusted"]="-n 0.02 -o 1 -e 5 -l 32"  
  ["permissive"]="-n 0.04 -o 2 -e 7 -l 100"  
)  

for prefix in "${!PARAMS[@]}"; do  
  singularity exec $CONTAINER bwa aln ${PARAMS[$prefix]} $GENOME test_R1.fq > ${prefix}_R1.sai  
  singularity exec $CONTAINER bwa sampe $GENOME ${prefix}_R1.sai ${prefix}_R2.sai test_R1.fq test_R2.fq > ${prefix}.sam  
done  
```

7. Final Recommendations
PoolSeq Priority: Favor specificity (low false positives) over sensitivity.

Validation:

Check allele frequencies at known loci.

Compare to GATK’s MarkDuplicates if PCR duplicates are a concern.

Citation:

bibtex
@article{kofler2016toolkit,
  title={Toolkit for PoolSeq analysis},
  author={Kofler, Robert and Schlötterer, Christian},
  journal={Genetics},
  year={2016}
}


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