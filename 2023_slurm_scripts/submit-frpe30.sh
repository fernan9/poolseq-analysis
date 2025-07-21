#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=SyncPileup-bulk    	#Name of the job


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

# pileup
singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
mpileup -f /scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Sorted/A1a_sorted.bam /scratch/frpe222/Sorted/A2a_sorted.bam \
/scratch/frpe222/Sorted/B1a_sorted.bam /scratch/frpe222/Sorted/B2a_sorted.bam \
/scratch/frpe222/Sorted/C1a_sorted.bam /scratch/frpe222/Sorted/C2a_sorted.bam \
/scratch/frpe222/Sorted/D1a_sorted.bam /scratch/frpe222/Sorted/D2a_sorted.bam  > /scratch/frpe222/PileupSync/PopAD.pileup

singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
java -ea -Xmx7g -jar /usr/local/popoolation2_1201/mpileup2sync.jar \
--input /scratch/frpe222/PileupSync/PopAD.pileup \
--output /scratch/frpe222/PileupSync/PopAD_java.sync \
--fastq-type sanger --min-qual 20 --threads 8