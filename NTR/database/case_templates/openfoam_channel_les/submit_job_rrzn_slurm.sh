#!/bin/bash -l
#SBATCH --job-name=<var CLUSTERJOBNAME var>
#SBATCH --nodes=<var NODES var>
#SBATCH --ntasks-per-node=<var TASKS_PER_NODE var>
#SBATCH --mem-per-cpu=2G
#SBATCH --time=<var WALLTIME var>
#SBATCH --constraint=[skylake|haswell]
#SBATCH --mail-user=<var JOB_EMAIL var>
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output batch.out
#SBATCH --error batch.err

module load GCC/4.9.3-2.25 OpenMPI/1.10.2 OpenFOAM/v1612+
source $FOAM_BASH
cd $SLURM_SUBMIT_DIR

mpirun  -np <var PROCS var> pimpleFoam -case ./ -parallel > log.log
