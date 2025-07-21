#!/bin/bash

#SBATCH -t 48:00:0                           	#Time for the job to run
#SBATCH --job-name=P1-6-VarianceSliding    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# MinCov 50

# Tajima's D
declare -A pool_sizes
pool_sizes=([P1]=39 [P2]=42 [P3]=43 [P4]=41 [P5]=48 [P6]=46)

for i in P1 P2 P3 P4 P5 P6; do
  singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
  --input /scratch/frpe222/Pileup/${i}.pileup --output /scratch/frpe222/VarianceSliding/${i}.d.10k \
  --measure D --window-size 10000 --step-size 10000 --fastq-type 'sanger' \
  --min-count 2 --min-coverage 5 --max-coverage 400 --min-qual 20 --pool-size "${pool_sizes[$i]}"

  singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
  --input /scratch/frpe222/VarianceSliding/${i}.d --output /scratch/frpe222/VarianceSliding/${i}.d.pdf \
  --ylab D --chromosomes "X 2L 2R 3L 3R 4"
done

# Pi diversity

for i in P1 P2 P3 P4 P5 P6; do
  singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
  --input /scratch/frpe222/Pileup/${i}.pileup --output /scratch/frpe222/VarianceSliding/${i}.pi.10k \
  --measure pi --window-size 10000 --step-size 10000 --fastq-type 'sanger' \
  --min-count 2 --min-coverage 5 --max-coverage 400 --min-qual 20 --pool-size "${pool_sizes[$i]}"

  singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
  --input /scratch/frpe222/VarianceSliding/${i}.pi --output /scratch/frpe222/VarianceSliding/${i}.pi.pdf \
  --ylab pi --chromosomes "X 2L 2R 3L 3R 4"
done

# Theta diversity

for i in P1 P2 P3 P4 P5 P6; do
  singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
  --input /scratch/frpe222/Pileup/${i}.pileup --output /scratch/frpe222/VarianceSliding/${i}.theta.10k \
  --measure theta --window-size 10000 --step-size 10000 --fastq-type 'sanger' \
  --min-count 2 --min-coverage 5 --max-coverage 400 --min-qual 20 --pool-size "${pool_sizes[$i]}"

  singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
  perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
  --input /scratch/frpe222/VarianceSliding/${i}.theta --output /scratch/frpe222/VarianceSliding/${i}.theta.pdf \
  --ylab theta --chromosomes "X 2L 2R 3L 3R 4"
done
