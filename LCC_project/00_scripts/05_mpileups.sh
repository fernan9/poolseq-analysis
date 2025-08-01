#!/bin/bash
#SBATCH --time=50:00:00  
#SBATCH --job-name=05_mpileups
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/05_mpileups_%j.err 
#SBATCH --output=000_logs/05_mpileups_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G  # Added memory request

module load ccs/singularity

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed"
CONTAINER="/share/singularity/images/ccs/ngsTools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/6_mpileups"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

# Get sorted list of BAM files with and without 'P' in the name
mapfile -t FILES_INBRED < <(find "$INPUT_DIR" -name "*.bam" | grep '/.*P.*\.bam$' | sort)
mapfile -t FILES_OUTBRED < <(find "$INPUT_DIR" -name "*.bam" | grep -v '/.*P.*\.bam$' | sort)

# Join arrays into space-separated strings
MPILEUP_FILE_INBRED="$OUTPUT_DIR/P16.mpileup"
MPILEUP_FILE_OUTBRED="$OUTPUT_DIR/ABCD.mpileup"

# To control sample names in samtools mpileup, the cleanest way is to use the -b option, which takes a file with a list of BAM file paths—one per line, in the order you want the samples to appear.
# However, samtools assigns sample names based on the @RG (read group) SM: field inside each BAM file.

echo "===== Starting MPileUp Inbred ===== at $(date)"

singularity run --app samtools113 "$CONTAINER" samtools \
    mpileup -B -q 20 -Q 20 -f "$GENOME" \
    "${FILES_INBRED[@]}" > "$MPILEUP_FILE_INBRED"

echo "===== Completed MPileUp Inbred ===== at $(date)"

echo "===== Starting MPileUp Outbred ===== at $(date)"

singularity run --app samtools113 "$CONTAINER" samtools \
    mpileup -B -q 20 -Q 20 -f "$GENOME" \
    "${FILES_OUTBRED[@]}" > "$MPILEUP_FILE_OUTBRED"

echo "===== Completed MPileUp Outbred ===== at $(date)"
