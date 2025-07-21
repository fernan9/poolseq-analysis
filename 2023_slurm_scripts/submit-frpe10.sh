#!/bin/bash

#SBATCH -t 48:00:0                           	#Time for the job to run
#SBATCH --job-name=A1a-VarianceSliding-MinCov    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node

#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# MinCov 50

# Tajima's D
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.d.15 \
--measure D --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 2 --min-coverage 15 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.d.15 --output /scratch/frpe222/VarianceSliding/A1a.d.15.pdf \
--ylab D --chromosomes "X 2L 2R 3L 3R 4"

# Pi diversity
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.pi.50 \
--measure pi --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 2 --min-coverage 50 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.pi.50 --output /scratch/frpe222/VarianceSliding/A1a.pi.50.pdf \
--ylab pi --chromosomes "X 2L 2R 3L 3R 4"

# Theta diversity
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.theta.50 \
--measure theta --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 2 --min-coverage 50 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.theta.50 --output /scratch/frpe222/VarianceSliding/A1a.theta.50.pdf \
--ylab theta --chromosomes "X 2L 2R 3L 3R 4"

# MinCov 100

# popoolation
# Tajima's D
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.d.10 \
--measure D --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 5 --min-coverage 10 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.d.10 --output /scratch/frpe222/VarianceSliding/A1a.d.10.pdf \
--ylab D --chromosomes "X 2L 2R 3L 3R 4"

# Pi diversity
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.pi.100 \
--measure pi --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 5 --min-coverage 100 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.pi.100 --output /scratch/frpe222/VarianceSliding/A1a.pi.100.pdf \
--ylab pi --chromosomes "X 2L 2R 3L 3R 4"

# Theta diversity
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.theta.100 \
--measure theta --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 5 --min-coverage 100 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.theta.100 --output /scratch/frpe222/VarianceSliding/A1a.theta.100.pdf \
--ylab theta --chromosomes "X 2L 2R 3L 3R 4"
