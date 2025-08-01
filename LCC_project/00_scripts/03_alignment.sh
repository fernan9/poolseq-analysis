#!/bin/bash
#SBATCH --time=30:00:00  
#SBATCH --job-name=03_alignment
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/03_alignment_%j.err 
#SBATCH --output=000_logs/03_alignment_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  
#SBATCH --cpus-per-task=8

module load ccs/singularity
# singularity run /share/singularity/images/ccs/popoolation/lcc-popoolation-ubuntu2204.sinf perl /opt/popoolation/basic-pipeline/identify-genomic-indel-regions.pl

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/3_trimmed"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/4_aligned"
CONTAINER="/share/singularity/images/ccs/conda/lcc-conda-2-centos8.sinf"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

# Get file list of all _1.fq files
FILES=($(ls ${INPUT_DIR}/*_R1.trimmed.fq.gz | sort))

echo "Starting script for alignment at $(date)"
echo "======================================"

for FILE1 in "${FILES[@]}"; do
    FILE2="${FILE1/_R1.trimmed.fq.gz/_R2.trimmed.fq.gz}"
    BASENAME=$(basename "$FILE1" | sed 's/_R1\.trimmed\.fq\.gz$//')
    ALIGNMENT_FILE1="${OUTPUT_DIR}/${BASENAME}_R1.sai"
    ALIGNMENT_FILE2="${OUTPUT_DIR}/${BASENAME}_R2.sai"
    OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}.sam"
    
    # check for paired file
    if [[ ! -f "$FILE2" ]]; then
        echo "Error: Paired file $FILE2 missing!" >&2
        continue
    fi
    # check for previous runs per SAMPLE
    if [[ -f "${ALIGNMENT_FILE1}" ]]; then
        echo "Skipping $BASENAME --- already aligned"
        continue
    else
        echo "Aligning $BASENAME --- missing alignment file"
        #continue
    fi

    echo "Aligning $BASENAME fwd"
    singularity run --app bwa0717 $CONTAINER bwa aln \
      -t $SLURM_CPUS_PER_TASK -o 1 -d 12 -e 5 -l 32 -n 0.02 \
      $GENOME $FILE1 > $ALIGNMENT_FILE1
    
    echo "Aligning $BASENAME rev"
    singularity run --app bwa0717 $CONTAINER bwa aln \
      -t $SLURM_CPUS_PER_TASK -o 1 -d 12 -e 5 -l 32 -n 0.02 \
      $GENOME $FILE2 > $ALIGNMENT_FILE2

    echo "Mapping $BASENAME"
    singularity run --app bwa0717 $CONTAINER bwa sampe \
      $GENOME $ALIGNMENT_FILE1 $ALIGNMENT_FILE2 \
      $FILE1 $FILE2 > $OUTPUT_FILE
    
    echo "Completed alignment for $BASENAME at $(date)"
    echo "Output: $OUTPUT_FILE"
    echo "======================================"
done

