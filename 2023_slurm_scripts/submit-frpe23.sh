#!/bin/bash

#SBATCH -t 76:00:0                           	#Time for the job to run
#SBATCH --job-name=Popool2-SyncPileup-P1-P6    	#Name of the job


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
for i in P1 P3 P4; do
  for j in P2 P5 P6; do
    # pileup
    singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
    mpileup -B /scratch/frpe222/Sorted/${i}_sorted.bam \
    /scratch/frpe222/Sorted/${j}_sorted.bam > /scratch/frpe222/PileupSync/${i}_${j}.pileup
  done
done

for i in P1 P3 P4; do
  for j in P2 P5 P6; do
    singularity run --app Popoolation2 /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
    java -ea -Xmx7g -jar /usr/local/popoolation2_1201/mpileup2sync.jar \
    --input /scratch/frpe222/PileupSync/${i}_${j}.pileup \
    --output /scratch/frpe222/PileupSync/${i}_${j}_java.sync \
    --fastq-type sanger --min-qual 20 --threads 8
  done
done