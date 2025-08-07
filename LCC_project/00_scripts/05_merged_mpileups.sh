#!/bin/bash
#SBATCH --time=50:00:00  
#SBATCH --job-name=05_merged_mpileups
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/05_merged_mpileups_%j.err 
#SBATCH --output=000_logs/05_merged_mpileups_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G  # Added memory request

module load ccs/singularity

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/5_processed"
CONTAINER="/share/singularity/images/ccs/ngsTools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
CONTAINER_SAMTOOLS="/share/singularity/images/ccs/conda/lcc-conda-10-rocky8.sinf"
CONTAINER_POP="/share/singularity/images/ccs/popoolation/lcc-popoolation-ubuntu2204.sinf"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/6_pileups"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

# --- Get merged BAMs (not SAMs!) ---
mapfile -t FILES < <(find "$INPUT_DIR" -name "merged*.bam" | sort)

echo "===== Starting Pipeline at $(date) ====="

for BAM_FILE in "${FILES[@]}"; do
    BASENAME=$(basename "$BAM_FILE" .bam)
    SORTED_BAM="$OUTPUT_DIR/${BASENAME}_sorted.bam"
    MPILEUP_FILE="$OUTPUT_DIR/${BASENAME}.mpileup"
    SYNC_FILE="$OUTPUT_DIR/${BASENAME}.sync"

    echo "Processing $BAM_FILE..."

    # --- Step 1: Sort Merged BAM ---
    if [[ -f "$SORTED_BAM" ]]; then
        echo "Skipping sorting bam, already exists at $SORTED_BAM."
    else
        singularity run --app samtools113 "$CONTAINER" \
            samtools sort -@ 8 -o "$SORTED_BAM" "$BAM_FILE" \
            || { echo "Sorting failed for $BAM_FILE"; exit 1; }

    # --- Step 2: Index Sorted BAM ---
        echo "Skipping sorting bam, already exists at $SORTED_BAM."

        singularity run --app samtools113 "$CONTAINER" \
            samtools index -@ 8 "$SORTED_BAM" \
            || { echo "Indexing failed for $SORTED_BAM"; exit 1; }
    fi

    # --- Step 3: Generate MPILEUP ---
    if [[ -f "$MPILEUP_FILE" ]]; then
        echo "Skipping mpilup, already exists at $MPILEUP_FILE."
    else
        singularity run --app samtools120 "$CONTAINER_SAMTOOLS" \
            samtools mpileup -d 1000 -B -f "$GENOME" \
            "$SORTED_BAM" -o $MPILEUP_FILE \
            || { echo "mpileup failed for $MPILEUP_FILE"; exit 1; }
    fi
    # --- Step 4: Convert to SYNC ---
    if [[ -f "$SYNC_FILE" ]]; then
        echo "Skipping sync, already exists at $SYNC_FILE."
    else
        singularity exec "$CONTAINER_POP" \
            java -ea -verbose:gc,class -jar /opt/popoolation2/mpileup2sync.jar \
            --input "$MPILEUP_FILE" \
            --output "$SYNC_FILE" \
            --fastq-type sanger \
            --min-qual 20 --threads 8 \
            || { echo "sync conversion failed for $MPILEUP_FILE"; exit 1; }
    fi
    echo "Completed $BASENAME at $(date)"
done

echo "===== Pipeline completed at $(date) ====="

singularity run /share/singularity/images/ccs/popoolation/lcc-popoolation-ubuntu2204.sinf perl /opt/popoolation2/fst-sliding.pl --help