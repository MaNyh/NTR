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

<var JOB_SOURCECMD var>

mpirun -mca btl '^openib' -np <var PROCS var> rhoPimpleFoam -parallel -case ./ > log.log
