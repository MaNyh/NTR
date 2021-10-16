#!/bin/bash -l
#SBATCH --job-name=<var JOB_NAME var>
#SBATCH --nodes=<var JOB_NODES var>
#SBATCH --ntasks-per-node=<var JOB_PPN var>
#SBATCH --mem-per-cpu=<var JOB_MEMPERCPU var>G
#SBATCH --time=<var JOB_RUNTIME var>
#SBATCH --mail-user=<var JOB_MAIL var>
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output batch.out
#SBATCH --error batch.err

<var JOB_SOURCECMD var>
cd $SLURM_SUBMIT_DIR

mpirun -np <var JOB_PROCS var> <var JOB_EXE var> -case ./ -parallel > log.log
