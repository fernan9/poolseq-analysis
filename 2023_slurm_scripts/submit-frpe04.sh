#!/bin/bash

#SBATCH -t 24:00:0                           	#Time for the job to run
#SBATCH --job-name=BWA-Mapping-Trimmed    	#Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user frpe222@uky.edu             #Where to send email

#SBATCH --account=coa_nmte222_uksr             	#Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#bwa 

# forward
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A0_trimed_1 > /scratch/frpe222/Trim/A0_trimed_1.sai

# reverse
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa aln -t 8 -o 2 -d 12 -e 12 -l 100 -n 0.01 \
/scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A0_trimed_2 > /scratch/frpe222/Trim/A0_trimed_2.sai

# mapping to sam
singularity run --app bwa0717 \
/share/singularity/images/ccs/conda/amd-conda1-centos8.sinf \
bwa sampe /scratch/frpe222/Genome/dmel-short-header.fa /scratch/frpe222/Trim/A0_trimed_1.sai /scratch/frpe222/Trim/A0_trimed_2.sai /scratch/frpe222/Trim/A0_trimed_1 /scratch/frpe222/Trim/A0_trimed_2 > /scratch/frpe222/Mapped/A0.sam

# "bwa aln" specifies that the "aln" algorithm within the BWA software package should be used for read alignment.
# "-t 8" specifies that the tool should use eight threads (i.e. eight CPU cores) for the alignment process, potentially speeding up the analysis if multiple cores are available.
# "-o 2" specifies the maximum number of gap opens allowed in the alignment process. In this case, up to two gap opens are permitted.
# "-d 12" specifies the maximum gap extension penalty that can be assigned in the alignment process.
# "-e 12" specifies the maximum edit distance allowed between the read and the reference genome, in terms of the number of mismatches and gaps.
# "-l 100" specifies the seed length used for the alignment process. In this case, a seed length of 100 base pairs is used.
# "-n 0.01" specifies the maximum allowed mismatch rate between the read and the reference genome. In this case, a maximum mismatch rate of 0.01 (or 1%) is allowed.