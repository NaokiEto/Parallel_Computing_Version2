#!/bin/bash

# Number of total processes          
NUMPROCESSES=$1

# Number of Child Processes

NUMCHILD=$(($NUMPROCESSES-1))

# Number of Threads
NUMTHREADS=$2

# File name of VTK produced
FILENAMEVTK=$3

# Number of times the program is to be ran
NUMREPETITIONS=$4

# Just parts of file names
MPI="MPI_NO_FILES_"
THREAD="_THREADS"
TRIAL="ChildProcess_TrialNumber"
ZERO="0"

# Make a directory of this particular amount of processes
mkdir $NUMCHILD$MPI$NUMTHREADS$THREAD
cd $NUMCHILD$MPI$NUMTHREADS$THREAD

# Do the program the amount of desired repeated times
for i in `seq 1 $NUMREPETITIONS`;
    do
        # Make the directory for the particular trial number
        mkdir $MPI$NUMCHILD$THREAD$NUMTHREADS$TRIAL$ZERO$i
        
        # Go back to the MPI_files_Real directory
        cd ../

        # Grab the vtk file and copy it into the trial number directory
        for f in *.vtk
            do
                cp -v $f $NUMCHILD$MPI$NUMTHREADS$THREAD/$MPI$NUMCHILD$THREAD$NUMTHREADS$TRIAL$ZERO$i
            done

        # Go to the trial number directory
        cd $NUMCHILD$MPI$NUMTHREADS$THREAD/$MPI$NUMCHILD$THREAD$NUMTHREADS$TRIAL$ZERO$i

        # Run the code 
        mpirun -np "$NUMPROCESSES" ./../../build/ApplyingVtkMarchingCubes "$NUMTHREADS" "$FILENAMEVTK"

        # Remove the original vtk file
        rm $f
        
        # Let's go back to the directory for the particular amount of processes and repeat
        cd ../
    done

