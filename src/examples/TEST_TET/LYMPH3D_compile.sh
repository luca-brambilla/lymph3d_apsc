#!/bin/bash
#PBS -S /bin/bash
#

# This is a bash script, thus every row is a bash command
# unless the first character of the row is a "#" :
# in this case the row is a bash comment.
# However, if a row starts with "#PBS" then is a Pbs/Torque directive
# When submitted to Pbs/Torque, this file will be executed by a bash shell
# on one node only, so the user must take care of propagating the needed instructions
# on every node of the requested pool (if multiple)
#-----------------SETTING THE REQUEST FOR THE HARDWARE ALLOCATION----#
#
# Put here the requests for the desired hardware:
# number of nodes (nodes)
# cores per node (ppn)
# walltime, i.e. maximum time of execution (HH:MM:SS)
# queue (name_of_the_queue)
#
# ACTUAL QUEUE LIMITS ARE REPORTED IN FILE /etc/motd
#
#PBS -l nodes=1:ppn=20,walltime=24:00:00 -q gigat
# Set the job name
#PBS -N LYMPH3D_PB
# Set the output file and merge it to the sterr
#PBS -o out-hostname-XyZ-N1x1-qsub.txt
#PBS -j oe
#PBS -e out-hostname-XyZ-N1x1.txt

#------------------SETTING THE ENVIRONMENT----------------------------#
#
# Here goes any environmental definition, such as module loading.
# BEWARE that any command will be executed on one node only,
# so take care of propagating any definition to every allocated node.
# Environmental definitions can be alternatively set in the .bashrc file
# export MY_ENVIRONMENTAL_VARIABLE=value
# Setting of environmental module system here
# module load MY_ENVIRONMENTAL_MODULE
# Start the job in the current directory (PBS starts in the home folder)
cd ${PBS_O_WORKDIR}
#
# Preliminaries
hostname
ulimit -s unlimited # make sure we can put big arrays on the stack
date

# Let's translate the file ${PBS_NODEFILE}, used by the
# PBS environment, into something that can be read by mpiexec.
sort ${PBS_NODEFILE} | uniq -c | awk '{ printf("%s\n", $2); }' > mpd.nodes
#-------------------RUN THE EXECUTABLE---------------------------------#
#
# YOUR LAUNCH COMMAND BELOW completed with I/O redirection if needed
# 1° example
#my_executable_command > output_file.txt 2>&1
# 2° example, with the MPI library
mpirun -machinefile mpd.nodes -n 1 -x OMP_NUM_THREADS=1 -npernode /u/gorlezza/LYMPH3D_elasticity_pb/Lymph3D 1>
lymph3d.out 2>&1