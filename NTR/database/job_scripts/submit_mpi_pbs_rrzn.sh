#PBS -S /bin/ksh
#PBS -N <var JOB_NAME var>
#PBS -j oe
#PBS -l nodes=<var JOB_NODES var>:ppn=<var JOB_PPN var>
#PBS -l walltime=<var JOB_RUNTIME var>,mem=16Gb
#PBS -M <var JOB_MAIL var>
#PBS -m ae
#PBS -W x=PARTITION:taurus:haku:lena
<var JOB_SOURCECMD var>
cd $PBS_O_WORKDIR
HOST=$(paste -s -d ',' ${PBS_NODEFILE})
mpirun --hostfile ${PBS_NODEFILE} -np <var JOB_PROCS var> <var JOB_EXE var> -parallel > log.log

