#!/bin/bash
#SBATCH --time=6:00:00  
#SBATCH --job-name=02_trimming
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/02_trimming_%j.err 
#SBATCH --output=000_logs/02_trimming_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  

module load ccs/singularity
# singularity run /share/singularity/images/ccs/popoolation/lcc-popoolation-ubuntu2204.sinf perl /opt/popoolation/basic-pipeline/identify-genomic-indel-regions.pl

# Directories
INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/1_raw_data"
OUTPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/3_trimmed"
CONTAINER="/share/singularity/images/ccs/popoolation/lcc-popoolation-ubuntu2204.sinf"
SCRIPT="/opt/popoolation/basic-pipeline/trim-fastq.pl"

# Get file list of all _1.fq files
FILES=($(ls ${INPUT_DIR}/*_1.fq | sort))

for FILE1 in "${FILES[@]}"; do
    FILE2="${FILE1/_1.fq/_2.fq}"
    BASENAME=$(basename "$FILE1" | sed 's/_L[0-9]*_1.fq//')
    OUTPUT1="${OUTPUT_DIR}/${BASENAME}_R1.trimmed.fq"

    if [[ -f "${OUTPUT1}.gz" ]]; then
        echo "Skipping $BASENAME --- already trimmed"
        continue
    else
        echo "Trimming $BASENAME --- missng file"
        #continue
    fi

    echo "Trimming $BASENAME"
    singularity run "$CONTAINER" perl "$SCRIPT" \
      --input1 "$FILE1" \
      --input2 "$FILE2" \
      --quality-threshold 20 \
      --min-length 50 \
      --fastq-type 'sanger' \
      --output1 "${OUTPUT_DIR}/${BASENAME}_R1.trimmed.fq" \
      --output2 "${OUTPUT_DIR}/${BASENAME}_R2.trimmed.fq" \
      --outputse "${OUTPUT_DIR}/${BASENAME}_SE.trimmed.fq"
done
