#!/bin/bash

#SBATCH -t 01:00:0                           	#Time for the job to run
#SBATCH --job-name=BWA-index-reference    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#bwa 

singularity run --app bwa0717 /share/singularity/images/ccs/conda/amd-conda1-centos8.sinf bwa index /scratch/frpe222/Genome/dmel-short-header.fa