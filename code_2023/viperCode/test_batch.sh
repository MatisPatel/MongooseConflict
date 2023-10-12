#!/bin/bash
#SBATCH -J mongoose_test
#SBATCH -n 1
#SBATCH -o logs/%N.%j.out
#SBATCH -e logs/%N.%j.err
#SBATCH -p compute
#SBATCH -t 00:20:00
#SBATCH --array=1-1000

# create NEW_TASK_ID which is the SLURM_ARRAY_TASK_ID but offset by a multiplier variable to cover all params


module purge
module add julia 

julia ../code/hpc_run.jl $SLURM_ARRAY_TASK_ID > logs/log.$SLURM_ARRAY_TASK_ID.txt

