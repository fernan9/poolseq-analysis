#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=Coverage   	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# popoolation2

for i in P1 P2 P3 P4 P5 P6; do
  # cross check
  singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
  coverage /scratch/frpe222/Sorted/${i}_sorted.bam -o /scratch/frpe222/Coverage/${i}.coverage
done

for i in A1a A2a B1a B2a C1a C2a D1a D2a; do
  # cross check
  singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
  coverage /scratch/frpe222/Sorted/${i}_sorted.bam -o /scratch/frpe222/Coverage/${i}.coverage
done