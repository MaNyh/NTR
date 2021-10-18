#!/bin/bash
###############################################################################
###                              job_TRACE.sh                               ###
###       start vom run-directory auf dem RRZN mit TRACE 9.3         	    ###
###                      10 July 2009 Florian Herbst                        ###
###                 modified by Hendrik Seehausen & Kenan Cengiz   May 2020 ###
###############################################################################
EMAIL=<var JOB_MAIL var>
NODES=<var JOB_NODES var>
NPROCS=<var JOB_PROCS var>
MEMORY=<var JOB_MEM var>
QUEUE="all"
WALLTIME=<var JOB_RUNTIME var>
#single or double precision TRACE version: "sp" or "dp"?
PRECISION="sp"
#Architecture is optional. Either leave blank or: haswell,nehalem,sandybridge,skylake
ARCHITECTURE=""
################################################################################
#DON'T CHANGE PART BELOW
################################################################################
BASENAME="TRACE.cgns"
#CONTROLFILE="TRACE_control.input"

let CPUS=($NPROCS \* $NODES)
BALANCE_FILE="BALANCE_"$CPUS"PROC"

cd input/
RUNDIR=$PWD
cd ../
WORKDIR=$PWD

#Jobname entspricht dem Ordnernamen
JOBNAME=${PWD##*/}

cd $RUNDIR

#append ":" before ARCHITECTURE if not empty
if [ "$ARCHITECTURE" != "" ]; then
   ARCHITECTURE=:$ARCHITECTURE
fi

#
#################################################################################
#Write tempscript to start TRACE
rm -f $RUNDIR/tempscript.sh
touch $RUNDIR/tempscript.sh
chmod +x $RUNDIR/tempscript.sh
echo -e "#!/bin/bash" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -N $JOBNAME" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -M $EMAIL" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -m abe" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -j oe" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -l nodes=$NODES$ARCHITECTURE:ppn=$NPROCS" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -l walltime=$WALLTIME" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -l mem=$MEMORY" >> $RUNDIR/tempscript.sh;
echo -e "#PBS -q $QUEUE" >> $RUNDIR/tempscript.sh;
echo -e "module purge" >> $RUNDIR/tempscript.sh;
echo -e "source $HOME/TRACE/trace94_${PRECISION}_profile.sh" >> $RUNDIR/tempscript.sh;

#Special MPI configuration for the enos cluster:
echo -e "[[ $HOSTNAME =~ ^enos-.* ]] && export MPI_IB_PKEY=0x8001" >> $RUNDIR/tempscript.sh;

#Using PREP
#echo -e "echo $CPUS" >> $RUNDIR/tempscript.sh;
echo -e "cd $WORKDIR/input/" >> $RUNDIR/tempscript.sh;
echo -e "PREP -cgns $BASENAME -clb -np $CPUS" >> $RUNDIR/tempscript.sh;

#Running trace
echo -e "cd $RUNDIR" >> $RUNDIR/tempscript.sh;
##### Mit Controle-File
echo -e "mpirun -d -machinefile \$PBS_NODEFILE TRACE -cgns $WORKDIR/input/$BASENAME -lb $WORKDIR/input/$BALANCE_FILE -o TRACE.lst." >> $RUNDIR/tempscript.sh;

#################################################################################
#submit job to queue
qsub $RUNDIR/tempscript.sh
