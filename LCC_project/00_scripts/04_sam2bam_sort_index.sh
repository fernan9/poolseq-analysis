#!/bin/bash
#SBATCH --time=50:00:00  
#SBATCH --job-name=04_sam2bam_sort_index
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/04_sam2bam_sort_index_%j.err 
#SBATCH --output=000_logs/04_sam2bam_sort_index_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G  # Added memory request

module load ccs/singularity

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/4_aligned"
CONTAINER="/share/singularity/images/ccs/ngsTools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed"

# Get sorted list of SAM files
mapfile -t FILES < <(find "$INPUT_DIR" -name "*.sam" | sort)

for SAM_FILE in "${FILES[@]}"; do
    REL_PATH="${SAM_FILE#$INPUT_DIR/}"  # removes input dir prefix
    BASENAME="${REL_PATH%.sam}"
    BAM_FILE="$OUTPUT_DIR/${BASENAME}.bam"

    echo "===== Starting $SAM_FILE ===== $(date)"

    # 1. Convert SAM to sorted BAM (save it for reuse)
    if [[ -f "$BAM_FILE" && -f "${BAM_FILE}.bai" ]]; then
        echo "Skipping BAM conversion, $BAM_FILE already exists."
    else        
        echo "Converting SAM to sorted BAM..."
        singularity run --app samtools113 "$CONTAINER" samtools view -@ $SLURM_CPUS_PER_TASK -b "$SAM_FILE" | \
        singularity run --app samtools113 "$CONTAINER" samtools sort -@ $SLURM_CPUS_PER_TASK -o "$BAM_FILE" - && \
        singularity run --app samtools113 "$CONTAINER" samtools index "$BAM_FILE"
    fi

    # 2. Run flagstat on the BAM file
    FLAGSTAT_FILE="$OUTPUT_DIR/${BASENAME}.flagstat.txt"
    if [[ -f "$FLAGSTAT_FILE" ]]; then
        echo "Skipping flagstat, already exists at $FLAGSTAT_FILE."
    else
        echo "Running flagstat..."
        singularity run --app samtools113 "$CONTAINER" samtools flagstat -@ $SLURM_CPUS_PER_TASK "$BAM_FILE" > "$FLAGSTAT_FILE"
    fi

    # 3. Calculate mean depth
    DEPTH_FILE="$OUTPUT_DIR/${BASENAME}.depth.txt"
    if [[ -f "$DEPTH_FILE" ]]; then
        echo "Skipping depth calculation, already exists at $DEPTH_FILE."
    else
        echo "Calculating depth metrics..."
        singularity run --app samtools113 "$CONTAINER" samtools depth -a "$BAM_FILE" | \
        awk '{sum+=$3; sumsq+=$3*$3} END {print "Mean depth:", sum/NR; print "Std dev:", sqrt(sumsq/NR - (sum/NR)^2)}' > "$DEPTH_FILE"
    fi
    
    echo "===== Completed $SAM_FILE ===== at $(date)"
done


