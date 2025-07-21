#!/bin/bash

#SBATCH -t 76:00:0                              #Time for the job to run
#SBATCH --job-name=Vcf_stats     #Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr              #Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "


# Define input and output directories
input_dir=/scratch/frpe222/Sorted
output_dir=/scratch/frpe222/Variant

for i in A1a A2a B1a B2a C1a C2a D1a D2a P1 P2 P3 P4 P5 P6; do
  # Run commands
  singularity run --app bcftools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf \
  bcftools stats "${output_dir}"/${i}_var.flt.vcf > "${output_dir}"/${i}_var.stats.txt
done