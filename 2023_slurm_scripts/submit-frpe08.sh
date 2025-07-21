#!/bin/bash

#SBATCH -t 08:00:0                           	#Time for the job to run
#SBATCH --job-name=Popoolation-variance-sliding    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#bwa 

# Tajima's D
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A0.pileup --output /scratch/frpe222/VarianceSliding/A0.d \
--measure D --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A0.d --output /scratch/frpe222/VarianceSliding/A0.d.pdf \
--ylab D --chromosomes "X 2L 2R 3L 3R 4"

# Pi diversity
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A0.pileup --output /scratch/frpe222/VarianceSliding/A0.pi \
--measure pi --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A0.pi --output /scratch/frpe222/VarianceSliding/A0.pi.pdf \
--ylab pi --chromosomes "X 2L 2R 3L 3R 4"