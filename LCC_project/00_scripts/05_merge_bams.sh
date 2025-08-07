#!/bin/bash
#SBATCH --time=50:00:00  
#SBATCH --job-name=05_merge_bams
#SBATCH --ntasks=1 
#SBATCH --partition=SKY32M192_L 
#SBATCH --error=000_logs/05_merge_bams_%j.err 
#SBATCH --output=000_logs/05_merge_bams_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G  # Added memory request

module load ccs/singularity

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed"
CONTAINER="/share/singularity/images/ccs/ngsTools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/6_pileups"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

# A1 - lib1 - header1.sam

echo "Merging A1 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header1.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_A1.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A0_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A1b_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A1a_H7MY2DSX5.bam" \
    || { echo "--- merging A1 set library 1 failed ---"; exit 1; }

# A2 - lib2 - header2.sam

echo "Merging A2 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header2.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_A2.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A0_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A2b_H7N25DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A2a_H7MY2DSX5.bam" \
    || { echo "--- merging A2 set library 1 failed ---"; exit 1; }

# B1 - lib3 - header3.sam

echo "Merging B1 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header3.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_B1.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B0_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B1b_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B1a_H7MY2DSX5.bam" \
    || { echo "--- merging B1 set library 1 failed ---"; exit 1; }

# B2 - lib4 - header4.sam

echo "Merging B2 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header4.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_B2.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B0_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B2b_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B2a_H7MY2DSX5.bam" \
    || { echo "--- merging B2 set library 1 failed ---"; exit 1; }

# C1 - lib5 - header5.sam

echo "Merging C1 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header5.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_C1.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C0_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C1b_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C1a_H7MY2DSX5.bam" \
    || { echo "--- merging C1 set library 1 failed ---"; exit 1; }

# C2 - lib6 - header6.sam

echo "Merging C2 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header6.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_C2.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C0_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C2b_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C2a_H7MY2DSX5.bam" \
    || { echo "--- merging C2 set library 1 failed ---"; exit 1; }

# D1 - lib7 - header7.sam

echo "Merging D1 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header7.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_D1.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D0_H7N25DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D1b_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D1a_H7MY2DSX5.bam" \
    || { echo "--- merging D1 set library 1 failed ---"; exit 1; }

# D2 - lib8 - header8.sam

echo "Merging D2 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header8.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_D2.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D0_H7N25DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D2b_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D2a_H7MY2DSX5.bam" \
    || { echo "--- merging D2 set library 1 failed ---"; exit 1; }

# ABCD1 - lib9 - header9.sam

echo "Merging ABCD1 set library 1"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header9.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_ABCD1.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A1a_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/A2a_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B1a_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/B2a_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C1a_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/C2a_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D1a_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/D2a_H7MY2DSX5.bam" \
    || { echo "--- merging ABCD1 set library 1 failed ---"; exit 1; }
