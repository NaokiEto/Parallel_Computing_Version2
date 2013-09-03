#!/bin/bash

# File name of VTK produced
FILENAMEVTK=$1

# File name of vtk to read
READFILE=$2

# Number of times the program is to be ran
NUMREPETITIONS=$3

# Just parts of file names
SERIAL="SERIAL"
TRIAL="_TrialNumber"
ZERO="0"

# Make a directory for testing Serial
mkdir $SERIAL
cd $SERIAL

# Do the program the amount of desired repeated times
for i in `seq 1 $NUMREPETITIONS`;
    do
        # Make the directory for the particular trial number
        mkdir $SERIAL$TRIAL$ZERO$i
        
        # Go back o the Serial_No_Files_Real directory
        cd ../

        # Grab the vtk file and copy it into the trial number directory
        for f in *.vtk
            do
                cp -v $f $SERIAL/$SERIAL$TRIAL$ZERO$i
            done

        # Go to the trial number directory
        cd $SERIAL/$SERIAL$TRIAL$ZERO$i

        # Run the code 
        ./../../build/ApplyingVtkMarchingCubes "$FILENAMEVTK" "$READFILE"

        # Remove the original vtk file
        rm $f
        
        # Let's go back to the directory for testing Serial and repeat
        cd ../
    done

