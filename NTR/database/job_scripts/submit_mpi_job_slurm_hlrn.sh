#!/bin/bash
#SBATCH -A <var JOB_ACCOUNT var>
#SBATCH --job-name=<var JOB_NAME var>
#SBATCH -t <var JOB_RUNTIME var>
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<var JOB_MAIL var>
#SBATCH --nodes=<var JOB_NODES var>
#SBATCH --tasks-per-node=<var JOB_PPN var>
#SBATCH -p standard96

export TMPDIR=$LOCAL_TMPDIR

<var JOB_SOURCECMD var>

mpirun -mca btl '^openib' -np <var JOB_PROCS var> <var JOB_EXE var> -parallel -case ./ > log.log
