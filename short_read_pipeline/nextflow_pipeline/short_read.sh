#!/bin/bash
#SBATCH --account=bgmp                      ### SLURM account which will be charged for the Job
#SBATCH --partition=bgmp                    ### Partition (like a queue in PBS)
#SBATCH --job-name=short_read                ### Job Name
#SBATCH --output=outfiles                   ### File in which to store job output
#SBATCH --error=errfiles                    ### File in which to store job error messages
#SBATCH --time=0-24:00:00                    ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                           ### Node count required for the job
#SBATCH --ntasks-per-node=1                 ### Nuber of tasks to be launched per Node
#SBATCH --cpus-per-task=16                   ### Number of cpus (cores) per task

#activate the nextflow environment
mamba activate /gpfs/projects/bgmp/shared/groups/2025/chimera/envs/Nextflow_env

#run nextflow script
/usr/bin/time -v nextflow run short_read.nf 