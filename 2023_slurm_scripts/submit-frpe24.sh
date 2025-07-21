#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=Popool2-Allelefreq-Fst-P1-6   	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# popoolation2

declare -A pool_sizes
pool_sizes=([P1_P2]="39:42" [P1_P5]="39:48" [P1_P6]="39:46" [P3_P2]="43:42" [P3_P5]="43:48" [P3_P6]="43:46" [P4_P2]="41:42" [P4_P5]="41:48" [P4_P6]="41:46")


for i in P1_P2 P1_P5 P1_P6 P3_P2 P3_P5 P3_P6 P4_P2 P4_P5 P4_P6; do
  singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation2_1201/snp-frequency-diff.pl \
  --input /scratch/frpe222/PileupSync/${i}_java.sync \
  --output-prefix /scratch/frpe222/PileupSync/${i} \
  --min-count 6 --min-coverage 50 --max-coverage 200
done



for i in P1_P2 P1_P5 P1_P6 P3_P2 P3_P5 P3_P6 P4_P2 P4_P5 P4_P6; do
  singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation2_1201/fst-sliding.pl \
  --input /scratch/frpe222/PileupSync/${i}_java.sync \
  --output /scratch/frpe222/PileupSync/${i}.fst \
  --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 \
  --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size "${pool_sizes[$i]}"
done


for i in P1_P2 P1_P5 P1_P6 P3_P2 P3_P5 P3_P6 P4_P2 P4_P5 P4_P6; do
  singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation2_1201//fst-sliding.pl \
  --input /scratch/frpe222/PileupSync/${i}_java.sync \
  --output /scratch/frpe222/PileupSync/${i}_w10000.fst \
  --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 \
  --min-covered-fraction 1 --window-size 10000 --step-size 10000 --pool-size "${pool_sizes[$i]}"
done