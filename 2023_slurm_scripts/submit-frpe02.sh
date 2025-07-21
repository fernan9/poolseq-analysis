#!/bin/bash

#SBATCH -t 01:00:0                           	#Time for the job to run
#SBATCH --job-name=Trimmomatic-test-PHRED20    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#Trimmomatic line execution

singularity run --app trimmomatic039 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf \
java -jar trimmomatic039 PE -trimlog ./trimmomatic-test-PHRED20.txt \
./RawData/A0_H7MY2DSX5_L1_1.fq.gz ./RawData/A0_H7MY2DSX5_L1_2.fq.gz \
./Trim/A0_forward_paired.fq.gz ./Trim/A0_forward_unpaired.fq.gz \
./Trim/A0_reverse_paired.fq.gz ./Trim/A0_reverse_unpaired.fq.gz \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40