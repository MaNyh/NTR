#PBS -S /bin/ksh
#PBS -N <var JOB_NAME var>
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=<var JOB_RUNTIME var>,mem=16Gb
#PBS -M <var JOB_MAIL var>
#PBS -m ae
#PBS -W x=PARTITION:taurus:haku:lena
<var JOB_SOURCECMD var>
cd $PBS_O_WORKDIR
<var JOB_EXE var> > log.log
