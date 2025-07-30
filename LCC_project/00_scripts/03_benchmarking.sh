#!/bin/bash
#SBATCH --time=12:00:00  
#SBATCH --job-name=03_benchmarking
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/03_benchmarking_%j.err 
#SBATCH --output=000_logs/03_benchmarking_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  
#SBATCH --cpus-per-task=8

module load ccs/singularity

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/3_trimmed"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/4_aligned"
CONTAINER="/share/singularity/images/ccs/conda/lcc-conda-2-centos8.sinf"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

# Test parameters
declare -A PARAMS=(
  ["strict"]="-n 0.01 -o 1 -e 3 -l 32"
  ["adjusted"]="-n 0.02 -o 1 -e 5 -l 32"
  ["permissive"]="-n 0.04 -o 2 -e 7 -l 100"
)

# Input files
FILE1="${INPUT_DIR}/B2a_H7MY2DSX5_R1.trimmed.fq.gz"
FILE2="${INPUT_DIR}/B2a_H7MY2DSX5_R2.trimmed.fq.gz"
BASENAME=$(basename "$FILE1" | sed 's/_R1\.trimmed\.fq\.gz$//')

# Verify input files exist
if [[ ! -f "$FILE1" || ! -f "$FILE2" ]]; then
  echo "Error: Input files missing! Check $FILE1 and $FILE2" >&2
  exit 1
fi

for prefix in "${!PARAMS[@]}"; do  
    ALIGNMENT_FILE1="${OUTPUT_DIR}/${prefix}_${BASENAME}_R1.sai"
    ALIGNMENT_FILE2="${OUTPUT_DIR}/${prefix}_${BASENAME}_R2.sai"
    OUTPUT_FILE="${OUTPUT_DIR}/${prefix}_${BASENAME}.sam"

    echo "===== Running $prefix parameters: ${PARAMS[$prefix]} ====="
    
    # Align forward reads
    echo "Aligning $BASENAME forward reads"
    singularity run --app bwa0717 "$CONTAINER" bwa aln \
      -t "$SLURM_CPUS_PER_TASK" -d 12 ${PARAMS[$prefix]} \
      "$GENOME" "$FILE1" > "$ALIGNMENT_FILE1" || { echo "bwa aln failed for $FILE1"; exit 1; }
    
    # Align reverse reads
    echo "Aligning $BASENAME reverse reads"
    singularity run --app bwa0717 "$CONTAINER" bwa aln \
      -t "$SLURM_CPUS_PER_TASK" -d 12 ${PARAMS[$prefix]} \
      "$GENOME" "$FILE2" > "$ALIGNMENT_FILE2" || { echo "bwa aln failed for $FILE2"; exit 1; }

    # Generate SAM
    echo "Generating SAM for $BASENAME"
    singularity run --app bwa0717 "$CONTAINER" bwa sampe \
      "$GENOME" "$ALIGNMENT_FILE1" "$ALIGNMENT_FILE2" \
      "$FILE1" "$FILE2" > "$OUTPUT_FILE" || { echo "bwa sampe failed"; exit 1; }
    
    echo "Completed $prefix alignment for $BASENAME at $(date)"
    echo "Output: $OUTPUT_FILE"
    echo "======================================"
done 
