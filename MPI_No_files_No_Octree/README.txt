This directory attempts to do MPI on marching cubes using the following
steps:

I find the vtk file, start MPI, divide up the vtk file, convert each piece to mhd, pass it through vtkMarchingCubes, send the resulting polydata to the master process (rank 0), which waits for all the processes. Then, the processes are conglomerated in serial and an vtk file is outputted.

To run this program, we would go 
to the PolyDataToImageData/build directory, add in the vtk file, and then
do the command

sudo mpirun -np 9 ./ApplyingVtkMarchingCubes AllStars.vtk

which includes the master process (so there are 8 child processes in this example)
