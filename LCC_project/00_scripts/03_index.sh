#!/bin/bash
#SBATCH --time=6:00:00  
#SBATCH --job-name=03_index
#SBATCH --ntasks=1 
#SBATCH --partition=CAL48M192_L 
#SBATCH --error=000_logs/03_index_%j.err 
#SBATCH --output=000_logs/03_index_%j.out 
#SBATCH --account=col_nmte222_uksr  
#SBATCH --mail-type ALL       
#SBATCH --mail-user fr_perezgalvez@outlook.com  

module load ccs/singularity
# singularity run /share/singularity/images/ccs/popoolation/lcc-popoolation-ubuntu2204.sinf perl /opt/popoolation/basic-pipeline/identify-genomic-indel-regions.pl

# Directories
CONTAINER="/share/singularity/images/ccs/conda/lcc-conda-2-centos8.sinf"
GENOME="/mnt/gpfs2_4m/scratch/frpe222/poolseq_project/0_genome/dmel_r6C.fasta"

singularity run --app bwa0717 $CONTAINER bwa index $GENOME