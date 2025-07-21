#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=FreqDiff-Fisher-Pop16   	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# popoolation2


singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation2_1201/snp-frequency-diff.pl \
--input /scratch/frpe222/PileupSync/Pop16_java.sync \
--output-prefix /scratch/frpe222/ExactFreq/P16 \
--min-count 2 --min-coverage 20 --max-coverage 2% 