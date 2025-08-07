#!/bin/bash
#SBATCH --time=50:00:00  
#SBATCH --job-name=05_merge_bam_p16
#SBATCH --ntasks=1 
#SBATCH --partition=SKY32M192_L 
#SBATCH --error=000_logs/05_merge_bam_p16_%j.err 
#SBATCH --output=000_logs/05_merge_bam_p16_%j.out 
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

echo "Merging P16 set library 10"
singularity run --app samtools113 "$CONTAINER" \
    samtools merge -r -h \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/header10.sam" -o \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/merged_P16.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/P1_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/P2_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/P3_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/P4_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/P5_H7MY2DSX5.bam" \
    "/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed/P6_H7MY2DSX5.bam" \
    || { echo "--- merging P16 set library 10 failed ---"; exit 1; }