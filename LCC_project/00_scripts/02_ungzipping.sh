#!/bin/bash
#SBATCH --time=24:00:00  
#SBATCH --job-name=02_ungzipping
#SBATCH --ntasks=1 
#SBATCH --partition=SKY32M192_L
#SBATCH --error=000_logs/02_ungzipping_%j.err 
#SBATCH --output=000_logs/02_ungzipping_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com

INPUT_DIR="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/1_raw_data"

cd "$INPUT_DIR" || exit

for f in *.fq.gz; do
  OUTFILE="${f%.gz}"
  if [ ! -f "$OUTFILE" ]; then
    echo "Decompressing $f"
    gunzip -c "$f" > "$OUTFILE"
  else
    echo "Skipping $f, $OUTFILE already exists"
  fi
done