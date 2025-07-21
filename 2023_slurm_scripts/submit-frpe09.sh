#!/bin/bash

#SBATCH -t 48:00:0                           	#Time for the job to run
#SBATCH --job-name=Sample-Pipeline-A1a    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

# trim fastqc popoolation
#singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf perl /usr/local/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
#--input1 /scratch/frpe222/RawData/A1a_H7MY2DSX5_L2_1.fq \
#--input2 /scratch/frpe222/RawData/A1a_H7MY2DSX5_L2_2.fq \
#--quality-threshold 20 --min-length 50 --fastq-type 'sanger' --output /scratch/frpe222/Trim/A1a_trimed

#bwa 

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A1a_trimed_1 > /scratch/frpe222/Trim/A1a_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A1a_trimed_2 > /scratch/frpe222/Trim/A1a_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A1a_trimed_1.sai /scratch/frpe222/Trim/A1a_trimed_2.sai \
/scratch/frpe222/Trim/A1a_trimed_1 /scratch/frpe222/Trim/A1a_trimed_2 > /scratch/frpe222/Mapped/A1a.sam

#filter sort bam

#samtools
singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
view -q 20 -b -S /scratch/frpe222/Mapped/A1a.sam > /scratch/frpe222/Mapped/A1a.bam
# sort
singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
sort /scratch/frpe222/Mapped/A1a.bam -o /scratch/frpe222/Sorted/A1a_sorted.bam

#samtools
# cross check
singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
flagstat /scratch/frpe222/Sorted/A1a_sorted.bam

# pileup
singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf samtools  \
mpileup /scratch/frpe222/Sorted/A1a_sorted.bam > /scratch/frpe222/Pileup/A1a.pileup

# popoolation
# Tajima's D
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.d \
--measure D --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.d --output /scratch/frpe222/VarianceSliding/A1a.d.pdf \
--ylab D --chromosomes "X 2L 2R 3L 3R 4"

# Pi diversity
singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Variance-sliding.pl \
--input /scratch/frpe222/Pileup/A1a.pileup --output /scratch/frpe222/VarianceSliding/A1a.pi \
--measure pi --window-size 50000 --step-size 10000 --fastq-type 'sanger' \
--min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 50

singularity run --app Popoolation /share/singularity/images/ccs/Popoolation/Popoolation.sinf \
perl /usr/local/popoolation_1.2.2/Visualise-output.pl \
--input /scratch/frpe222/VarianceSliding/A1a.pi --output /scratch/frpe222/VarianceSliding/A1a.pi.pdf \
--ylab pi --chromosomes "X 2L 2R 3L 3R 4"