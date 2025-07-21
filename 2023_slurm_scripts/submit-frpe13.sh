#!/bin/bash

#SBATCH -t 48:00:0                           	#Time for the job to run
#SBATCH --job-name=BWA-Mapping-Trimmed-P1-6    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#P1 

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P1_trimed_1 > /scratch/frpe222/Trim/P1_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P1_trimed_2 > /scratch/frpe222/Trim/P1_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P1_trimed_1.sai /scratch/frpe222/Trim/P1_trimed_2.sai \
/scratch/frpe222/Trim/P1_trimed_1 /scratch/frpe222/Trim/P1_trimed_2 > /scratch/frpe222/Mapped/P1.sam

#P2

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P2_trimed_1 > /scratch/frpe222/Trim/P2_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P2_trimed_2 > /scratch/frpe222/Trim/P2_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P2_trimed_1.sai /scratch/frpe222/Trim/P2_trimed_2.sai \
/scratch/frpe222/Trim/P2_trimed_1 /scratch/frpe222/Trim/P2_trimed_2 > /scratch/frpe222/Mapped/P2.sam

#P3

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P3_trimed_1 > /scratch/frpe222/Trim/P3_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P3_trimed_2 > /scratch/frpe222/Trim/P3_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P3_trimed_1.sai /scratch/frpe222/Trim/P3_trimed_2.sai \
/scratch/frpe222/Trim/P3_trimed_1 /scratch/frpe222/Trim/P3_trimed_2 > /scratch/frpe222/Mapped/P3.sam

#P4 

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P4_trimed_1 > /scratch/frpe222/Trim/P4_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P4_trimed_2 > /scratch/frpe222/Trim/P4_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P4_trimed_1.sai /scratch/frpe222/Trim/P4_trimed_2.sai \
/scratch/frpe222/Trim/P4_trimed_1 /scratch/frpe222/Trim/P4_trimed_2 > /scratch/frpe222/Mapped/P4.sam

#P5 

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P5_trimed_1 > /scratch/frpe222/Trim/P5_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P5_trimed_2 > /scratch/frpe222/Trim/P5_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P5_trimed_1.sai /scratch/frpe222/Trim/P5_trimed_2.sai \
/scratch/frpe222/Trim/P5_trimed_1 /scratch/frpe222/Trim/P5_trimed_2 > /scratch/frpe222/Mapped/P5.sam

#P6 

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P6_trimed_1 > /scratch/frpe222/Trim/P6_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P6_trimed_2 > /scratch/frpe222/Trim/P6_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa \
/scratch/frpe222/Trim/P6_trimed_1.sai /scratch/frpe222/Trim/P6_trimed_2.sai \
/scratch/frpe222/Trim/P6_trimed_1 /scratch/frpe222/Trim/P6_trimed_2 > /scratch/frpe222/Mapped/P6.sam