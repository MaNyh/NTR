#!/bin/bash -l
#SBATCH --job-name=<var CLUSTERJOBNAME var>
#SBATCH --nodes=<var NODES var>
#SBATCH --ntasks-per-node=<var TASKS_PER_NODE var>
#SBATCH --mem-per-cpu=2G
#SBATCH --time=100:00:00
#SBATCH --constraint=[skylake|haswell]
#SBATCH --mail-user=<var JOB_EMAIL var>
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output batch.out
#SBATCH --error batch.err

source ~/ModuleLoading/OpenFoam1612.sh
source $FOAM_BASH
cd $SLURM_SUBMIT_DIR

mpirun -np <var PROCS var> rhoPimpleFoam -case ./ -parallel > log.log
