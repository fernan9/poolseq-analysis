#!/bin/bash

#SBATCH -t 48:00:0                           	#Time for the job to run
#SBATCH --job-name=FastQC-multiple                	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#FastQC program execution command
for i in  /home/frpe222/RawData/*.fq.gz
do
	singularity run --app fastqc0119 /share/singularity/images/ccs/conda/amd-conda6-rocky8.sinf fastqc --output=/home/frpe222/RawData/Example-out $i
done