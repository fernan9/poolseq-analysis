#!/bin/bash
#SBATCH --time=01:00:00  
#SBATCH --job-name=02_ungzipping
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_D 
#SBATCH --error=000_logs/02_ungzipping_%j.err 
#SBATCH --output=000_logs/02_ungzipping_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com

INPUT_DIR="/project/nmte222_uksr/perezgalvez/poolseq_project/1_raw_data"


cd "$INPUT_DIR"

for f in *.fq.gz; do
  echo "Decompressing $f"
  gunzip -c "$f" > "${f%.gz}"
done