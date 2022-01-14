#!/bin/bash

#SBATCH --qos=medium
#SBATCH --job-name=Basin_Stability_DAE
#SBATCH --account=coen
#SBATCH --error=%x-%j-%N.err
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --tasks-per-node=10

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.6.1
julia ieee96.jl $SLURM_NTASKS
