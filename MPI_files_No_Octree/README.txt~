This directory attempts to do MPI on marching cubes using the following
steps:

get the VTK file,  divide up it up and do MPI, convert each piece to metaimage data, apply marching cubes to each process, output the resulting vtk data, and then conglomerate the files into 1 vtk file.

To run this program, we would go 
to the PolyDataToImageData/build directory, add in the vtk file, and then
do the command

sudo mpirun -np 9 ./ApplyingVtkMarchingCubes AllStars.vtk

which includes the master process (so there are 8 child processes in this example)



