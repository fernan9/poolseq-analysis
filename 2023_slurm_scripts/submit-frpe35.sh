#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=FST-PopAD				   	#Name of the job


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
perl /usr/local/popoolation2_1201/fst-sliding.pl \
--input /scratch/frpe222/PileupSync/PopAD_java.sync \
--output /scratch/frpe222/Fst/PAD.fst \
--min-count 2 --min-coverage 20 --max-coverage 2% --pool-size 23:40:43:34:21:35:46:46 \
--window-size 10000 --step-size 10000

echo "Job $SLURM_JOB_ID running (2nd PART) on SLURM NODELIST: $SLURM_NODELIST "

singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation2_1201/fst-sliding.pl \
--input /scratch/frpe222/PileupSync/PopAD_java.sync \
--output /scratch/frpe222/Fst/PAD_karlsson.fst \
--min-count 2 --min-coverage 20 --max-coverage 2% --pool-size 23:40:43:34:21:35:46:46 \
--window-size 10000 --step-size 10000 --karlsson-fst