#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=Popool2-Allelefreq-Fst    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# popoolation2

declare -A pool_sizes
pool_sizes=([A]=23 [B]=34 [C]=21 [D]=46)


for i in A B C D; do
  singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation2_1201/snp-frequency-diff.pl \
  --input /scratch/frpe222/PileupSync/${i}1a_${i}2a_java.sync \
  --output-prefix /scratch/frpe222/PileupSync/${i}1a_${i}2a \
  --min-count 6 --min-coverage 50 --max-coverage 200
done



for i in A B C D; do
  singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation2_1201//fst-sliding.pl \
  --input /scratch/frpe222/PileupSync/${i}1a_${i}2a_java.sync \
  --output /scratch/frpe222/PileupSync/${i}1a_${i}2a.fst \
  --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 \
  --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size "${pool_sizes[$i]}"
done
