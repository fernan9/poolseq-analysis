#!/bin/bash
#SBATCH --time=50:00:00  
#SBATCH --job-name=05_Ne_mpileups
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/05_Ne_mpileups_%j.err 
#SBATCH --output=000_logs/05_Ne_mpileupss_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G  # Added memory request

module load ccs/singularity

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed"
CONTAINER="/share/singularity/images/ccs/ngsTools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
CONTAINER_POP="/share/singularity/images/ccs/popoolation/lcc-popoolation-ubuntu2204.sinf"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/6_pileups"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

# Sample File names and location

FILE_NAMES_A1="$OUTPUT_DIR/names_A1.txt"

# Output names of mpileups

MPILEUP_A1="$OUTPUT_DIR/A1.mpileup"
SYNC_A1="$OUTPUT_DIR/A1.sync"

# To control sample names in samtools mpileup, the cleanest way is to use the -b option, which takes a file with a list of BAM file paths—one per line, in the order you want the samples to appear.
# However, samtools assigns sample names based on the @RG (read group) SM: field inside each BAM file.

echo "===== Starting MPileUp A1 at $(date) ====="

singularity run --app samtools113 "$CONTAINER" samtools \
    mpileup -B -q 20 -Q 20 -f "$GENOME" --threads 8 \
    -b "$FILE_NAMES_A1" -o "$MPILEUP_A1" \
    || { echo "mpileup failed"; exit 1; }

# Step 2: Convert to SYNC
singularity run --app Popoolation2 $CONTAINER_POP \
    perl /opt/popoolation2/mpileup2sync.pl \
    --input "$MPILEUP_A1" \
    --output "$SYNC_A1" \
    --fastq-type sanger --min-qual 20\
    || { echo "sync conversion failed"; exit 1; }

echo "===== Pipeline completed at $(date) ====="

singularity run --app samtools113 /share/singularity/images/ccs/ngsTools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools view -H A0_H7MY2DSX5.bam 
