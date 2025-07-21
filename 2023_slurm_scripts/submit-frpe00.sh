#!/bin/bash

#SBATCH -t 00:15:0                           	#Time for the job to run
#SBATCH --job-name=FastQC-test                	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#FastQC program execution command
singularity run --app fastqc0119 /share/singularity/images/ccs/conda/amd-conda6-rocky8.sinf fastqc A0_H7MY2DSX5_L1_1.fq.gz A0_H7MY2DSX5_L1_2.fq.gz A1a_H7MY2DSX5_L2_1.fq.gz A1a_H7MY2DSX5_L2_2.fq.gz