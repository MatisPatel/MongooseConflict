#!/bin/bash
#SBATCH -J mongoose_test
#SBATCH -n 1
#SBATCH -o logs/%N.%j.out
#SBATCH -e logs/%N.%j.err
#SBATCH -p compute
#SBATCH -t 00:20:00
#SBATCH --array=1-729

module purge
module add julia 

julia ../code/hpc_run.jl $SLURM_ARRAY_TASK_ID > logs/log.$SLURM_ARRAY_TASK_ID.txt

