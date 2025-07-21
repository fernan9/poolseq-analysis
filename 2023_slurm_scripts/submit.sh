#!/bin/bash

#SBATCH -t 00:05:0                              #Time for the job to run
#SBATCH --job-name=Conda-example                #Name of the job


#SBATCH -n 1                                    #Number of cores needed for the job; cannot be more than the number of cores available on a single node
#SBATCH --partition=normal                      #Name of the queue

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user gazula@uky.edu              #Where to send email

#SBATCH --account=coa_vgazu2_uksr               #Name of account to run under

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

#Velvet program execution command
singularity run --app velvet1210 /share/singularity/images/ccs/conda/lcc-bio-1-centos7.sinf velvetg

