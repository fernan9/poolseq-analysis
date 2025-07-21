#!/bin/bash
#SBATCH --time=08:15:00                 # Time limit for the job (REQUIRED).
#SBATCH --job-name=01_fastqc            # Job name
#SBATCH --ntasks=1                      # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --partition=SKY32M192_L         # Partition/queue to run the job in. (REQUIRED)
#SBATCH --error=000_logs/01_fastqc_%j.err   # Error file for this job.
#SBATCH --output=000_logs/01_fastqc_%j.out  # Output file for this job.
#SBATCH --account=col_nmte222_uksr      # Project allocation account name (REQUIRED)
#SBATCH --mail-type ALL                 # Send email on start/end
#SBATCH --mail-user fr_perezgalvez@outlook.com     # Where to send email

#Module needed to run the conda job
module load ccs/singularity

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#FastQC program execution command
for file in 1_raw_data/*.fq.gz; do
	singularity run --app fastqc0119 /share/singularity/images/ccs/conda/lcc-conda-2-centos8.sinf fastqc --output=/project/nmte222_uksr/perezgalvez/poolseq_project/2_fastqc/ $file
    echo "File $file done"
done
