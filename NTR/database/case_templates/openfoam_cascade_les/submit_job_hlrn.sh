#!/bin/bash
#SBATCH -A <var HLRN_JOB_ACCOUNT var>
#SBATCH --job-name=<var CLUSTERJOBNAME var>
#SBATCH -t <var RUNTIME var>
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<var JOB_EMAIL var>
#SBATCH --nodes=<var NODES var>
#SBATCH --tasks-per-node=<var TASKS_PER_NODE var>
#SBATCH -p standard96

export TMPDIR=$LOCAL_TMPDIR

module load  gcc/7.5.0 cmake/3.16.2; source /home/niikceng/OpenFOAM/OpenFOAM-v1612+/etc/bashrc
export FOAM_USER_PATH='$HOME/OpenFOAM/niinyhma-v1612+/'
mpirun -mca btl '^openib' -np <var PROCS var> rhoPimpleFoam -parallel -case ./ > log.log

