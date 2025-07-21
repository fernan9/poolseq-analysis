#!/bin/bash

#SBATCH -t 08:00:0                           	#Time for the job to run
#SBATCH --job-name=Crosscheck-Pileup   	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#samtools
# cross check
singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
flagstat /scratch/frpe222/Sorted/A0_sorted.bam

# pileup
singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
mpileup /scratch/frpe222/Sorted/A0_sorted.bam > /scratch/frpe222/Pileup/A0.pileup