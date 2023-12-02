#!/bin/bash
#SBATCH -J mongoose_test
#SBATCH -n 1
#SBATCH -o logs/%N.%j.out
#SBATCH -e logs/%N.%j.err
#SBATCH -p compute
#SBATCH -t 00:59:59
#SBATCH --array=1-1000

mult=$1
# create NEW_TASK_ID which is the SLURM_ARRAY_TASK_ID but offset by a multiplier variable to cover all params
NEW_TASK_ID=$( echo ${SLURM_ARRAY_TASK_ID} | awk '{ print $1+('$mult'*1000)}' )

module purge
module add julia 

julia ../code/hpc_run.jl $NEW_TASK_ID > logs/log.$NEW_TASK_ID.txt

