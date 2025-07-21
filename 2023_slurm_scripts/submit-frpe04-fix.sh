#!/bin/bash

#SBATCH -t 24:00:0                           	#Time for the job to run
#SBATCH --job-name=BWA-Mapping-to-SAM    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#bwa 

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A0_trimed_1.sai /scratch/frpe222/Trim/A0_trimed_2.sai /scratch/frpe222/Trim/A0_trimed_1 /scratch/frpe222/Trim/A0_trimed_2 > /scratch/frpe222/Mapped/A0.sam