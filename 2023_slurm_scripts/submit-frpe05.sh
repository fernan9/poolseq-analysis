#!/bin/bash

#SBATCH -t 08:00:0                           	#Time for the job to run
#SBATCH --job-name=Popoolation-trim    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#bwa 

# 
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl /usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/A0_H7MY2DSX5_L1_1.fq \
--input2 /scratch/frpe222/RawData/A0_H7MY2DSX5_L1_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/A0_trimed

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl /usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/A1a_H7MY2DSX5_L2_1.fq \
--input2 /scratch/frpe222/RawData/A1a_H7MY2DSX5_L2_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/A1a_trimed