#!/bin/bash

#SBATCH -t 08:00:0                           	#Time for the job to run
#SBATCH --job-name=Popoolation-trim-P1-6    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#bwa 

# P1
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl \
/usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/P1_H7MY2DSX5_L3_1.fq \
--input2 /scratch/frpe222/RawData/P1_H7MY2DSX5_L3_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/P1_trimed

# P2
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl \
/usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/P2_H7MY2DSX5_L3_1.fq \
--input2 /scratch/frpe222/RawData/P2_H7MY2DSX5_L3_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/P2_trimed

# P3
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl \
/usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/P3_H7MY2DSX5_L3_1.fq \
--input2 /scratch/frpe222/RawData/P3_H7MY2DSX5_L3_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/P3_trimed

# P4
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl \
/usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/P4_H7MY2DSX5_L3_1.fq \
--input2 /scratch/frpe222/RawData/P4_H7MY2DSX5_L3_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/P4_trimed

# P5
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl \
/usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/P5_H7MY2DSX5_L3_1.fq \
--input2 /scratch/frpe222/RawData/P5_H7MY2DSX5_L3_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/P5_trimed

# P6
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl \
/usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /scratch/frpe222/RawData/P6_H7MY2DSX5_L3_1.fq \
--input2 /scratch/frpe222/RawData/P6_H7MY2DSX5_L3_2.fq \
--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/P6_trimed