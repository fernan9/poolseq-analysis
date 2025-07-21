#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=Fish-Pop16   	#Name of the job


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

export PERL5LIB=/mnt/gpfs3_amd/share/apps/spack/0.16.2-3602-af806a8c1e/opt/spack/linux-centos8-zen2/gcc-9.3.0/perl-5.34.0-2kus2s4rxtdvqmzg6dctmgsztl7pnmc2/lib/site_perl/5.34.0

singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation2_1201/fisher-test.pl \
--input /scratch/frpe222/PileupSync/Pop16_java.sync \
--output /scratch/frpe222/Fisher/Pop16.fet \
--min-count 2 --min-coverage 20 --max-coverage 2%

echo "Job $SLURM_JOB_ID running (2nd PART) on SLURM NODELIST: $SLURM_NODELIST "

singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation2_1201/fisher-test.pl \
--input /scratch/frpe222/PileupSync/Pop16_java.sync \
--output /scratch/frpe222/Fisher/Pop16_w1000_geom.fet \
--min-count 2 --min-coverage 20 --max-coverage 2% --min-covered-fraction 0.2 \
--window-size 10000 --step-size 10000 --window-summary-method geometricmean

echo "Job $SLURM_JOB_ID running (2nd PART) on SLURM NODELIST: $SLURM_NODELIST "

singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation2_1201/fisher-test.pl \
--input /scratch/frpe222/PileupSync/Pop16_java.sync \
--output /scratch/frpe222/Fisher/Pop16_w1000_mult.fet \
--min-count 2 --min-coverage 20 --max-coverage 2% --min-covered-fraction 0.2 \
--window-size 10000 --step-size 10000