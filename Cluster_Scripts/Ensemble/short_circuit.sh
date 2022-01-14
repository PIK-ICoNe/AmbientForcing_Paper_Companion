#!/bin/bash

#SBATCH --qos=standby
#SBATCH --partition=priority
#SBATCH --account=coen
#SBATCH --error=%x-%j-%N.err
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --tasks-per-node=16
#SBATCH --job-name=short_circuit

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.6.1
julia short_circuit.jl $SLURM_NTASKS
