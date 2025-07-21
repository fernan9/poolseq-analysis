#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=Popool2-Fisher-AD    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node


#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# popoolation2

for i in A B C D; do
  singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation2_1201/fisher-test.pl \
  --input /scratch/frpe222/PileupSync/${i}1a_${i}2a_java.sync \
  --output /scratch/frpe222/PileupSync/${i}1a_${i}2a.fet \
  --min-count 6 --min-coverage 50 --max-coverage 200 --suppress-noninformative
done
