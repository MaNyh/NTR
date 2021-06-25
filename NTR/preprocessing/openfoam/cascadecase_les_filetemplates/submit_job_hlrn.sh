#!/bin/bash
#SBATCH -A niinyhma
#SBATCH --job-name=P14_start
#SBATCH -t __RUNTIME__
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nyhuis@tfd.uni-hannover.de
#SBATCH --nodes=__NODES__
#SBATCH --tasks-per-node=96
#SBATCH -p standard96

export TMPDIR=$LOCAL_TMPDIR

#source /home/niinyhma/OpenFOAM/OpenFOAM-v1612+/etc/bashrc

module load  gcc/7.5.0 cmake/3.16.2; source /home/niikceng/OpenFOAM/OpenFOAM-v1612+/etc/bashrc
export FOAM_USER_PATH='$HOME/OpenFOAM/niinyhma-v1612+/'
mpirun -mca btl '^openib' -np __PROCS__ rhoPimpleFoam -parallel -case ./ > log.log

