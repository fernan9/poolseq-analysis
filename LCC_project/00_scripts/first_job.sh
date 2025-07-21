#!/bin/bash
#SBATCH --time=00:15:00                 # Time limit for the job (REQUIRED).
#SBATCH --job-name=my_test_job          # Job name
#SBATCH --ntasks=1                      # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --partition=SKY32M192_D         # Partition/queue to run the job in. (REQUIRED)
#SBATCH --error=logs/slurm-%j.err                 # Error file for this job.
#SBATCH --output=logs/slurm-%j.out                 # Output file for this job.
#SBATCH --account=col_nmte222_uksr      # Project allocation account name (REQUIRED)
#SBATCH --mail-type ALL                 # Send email on start/end
#SBATCH --mail-user frpe222@uky.edu     # Where to send email

echo "Hello world. This is my first job"   # This is the program that will be executed. You will substitute this with your scientific program.