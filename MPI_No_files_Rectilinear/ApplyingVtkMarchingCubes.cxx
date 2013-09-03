/**
* Do whatever you want with public license
* Version 1, August 11, 2013
*
* Copyright (C) 2013 Naoki Eto <neto@lbl.gov>
*
* Everyone is permitted to copy and distribute verbatim or modified
* copies of this license document, and changing it is allowed as long
* as the name is changed.
*
* Do whatever you want with the public license
*
* TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
*
* 0. You just do what you want to do.
* 
*/
/**
* @file ApplyingVtkMarchingCubes.cxx
* @author Naoki Eto
* @date August 12, 2013
* @brief This program gets the VTK file,  divides up it up, and performs MPI. 
*        It convert each piece to metaimage data so that vtkMarchingCubes 
*        class can be applied, apply marching cubes to each process, outputs 
*        the resulting vtk data, and then conglomerates the data into 1 vtk 
*        file in the master process.
* @param[in] number of processes - number of processes for MPI (look at README 
             for more information)
* @param[in] argv[1] - the output's filename
* @param[out] pWriter - vtkPolyData file with the output's filename
* @return - EXIT_SUCCESS at the end
*/

#include <dirent.h>
#include <vector>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <string.h>
#include <vtkPolyDataMapper.h>
#include <vtkExtractVOI.h>
#include <vtkMarchingCubes.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDataArray.h>
#include <vtkLookupTable.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkMPIController.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>

#include <mpi.h>//"/home/users/neto/mpich2-1.4.1p1-install/include/mpi.h"
#include <stdio.h>
#include "/work2/vt-system-install/include/vampirtrace/vt_user.h"
#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPointSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkDataReader.h>
#include <vtkUnstructuredGrid.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkRectilinearGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGridGeometryFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkPoints.h>
#include <algorithm> 

#include <time.h>

void process(int procRank, int procSize, vtkMPIController* procController, const char* fp)
{
    /* The vtk file extension we want to search for */

    //VT_ON();
    //VT_USER_START("Region 1");

    int NumOfCharPD = strlen(fp) + 5 + 1;
    int original = NumOfCharPD;

    // figure out how many characters are in the temporary file name, 
    if(procRank - 1 > 9)
        NumOfCharPD += (int) log10(procRank - 1);

    char buf[(int) NumOfCharPD - original + 1];
 
    sprintf(buf, "%d", procRank-1);

    const char* suffix = ".vtk";

    const char* prefix = fp;

    char prefix_suffix[NumOfCharPD];

    strcpy(prefix_suffix, prefix);

    strcat(prefix_suffix, buf);

    strcat(prefix_suffix, suffix);

    //VT_USER_END("Region 1");
    //VT_OFF();


    vtkRectilinearGridReader *reader = vtkRectilinearGridReader::New();

    reader->SetFileName(prefix_suffix);
    reader->Update();

    // Create a grid
    vtkSmartPointer<vtkRectilinearGrid> grid = reader->GetOutput();

    vtkDataArray* dim;

    vtkPointData* pointdata = grid->GetPointData();

    double* range;

    range = pointdata->GetArray("grad")->GetRange();

    vtkContourFilter* contour = vtkContourFilter::New();

    // name of array is "grad"
    contour->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "grad");

    // better than setinput
    contour->SetInputConnection(grid->GetProducerPort());

    // woo 50 contours
    contour->GenerateValues(50, range);

    contour->ComputeNormalsOn();

    contour->Update();

    vtkPolyData *smoothed_polys = contour->GetOutput();

    // calc cell normal
    vtkPolyDataNormals *triangleCellNormals= vtkPolyDataNormals::New();

    #if VTK_MAJOR_VERSION <= 5
        triangleCellNormals->SetInputConnection(contour->GetOutputPort());
    #else
        triangleCellNormals->SetInputData(smoothed_polys);   
    #endif

    triangleCellNormals->ComputeCellNormalsOn();
    triangleCellNormals->ComputePointNormalsOff();
    triangleCellNormals->ConsistencyOn();
    triangleCellNormals->AutoOrientNormalsOn();
    triangleCellNormals->Update(); // creates vtkPolyData

    // send the vtkPolyData to the master process
    #if VTK_MAJOR_VERSION <= 5
        procController->Send(triangleCellNormals->GetOutput(), 0, 101);
    #else
        procController->Send(triangleCellNormals->GetOutputPort(), 0, 101);  
    #endif

/*

    //rendering functions

    vtkPolyDataMapper* contourMapper = vtkPolyDataMapper::New();
    contourMapper->SetInput(triangleCellNormals->GetOutput());

     vtkActor* contourActor = vtkActor::New();
     contourActor->SetMapper(contourMapper);

    // Create a renderer, render window, and interactor
    vtkRenderer* renderer = vtkRenderer::New();
    vtkRenderWindow *renderWindow = vtkRenderWindow::New();
    vtkRenderWindowInteractor *renderWindowInteractor = vtkRenderWindowInteractor::New();

    // Add the actors to the scene
    renderWindow->AddRenderer(renderer);
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderer->AddActor(contourActor);

    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();
*/
}

/**
 * This program converts a vtkPolyData image into volume representation 
 * (vtkImageData) where the foreground voxels are 1 and the background 
 * voxels are 0. Internally vtkPolyDataToImageStencil is utilized as 
 * as MPI. The resultant image is saved as metaimage data. vtkMarchingCubes 
 * is applied to these pieces of data, which are then conglomerated by the 
 * master process, and an output vtk file is outputted.
 */
int main(int argc, char *argv[])
{
    vtkMPIController* controller = vtkMPIController::New();

    double t1 = MPI_Wtime();

    // Initializing MPI
    controller->Initialize(&argc, &argv);

    /* Figure out total amount of processors */
    int rank = controller->GetLocalProcessId();

    /* Figure out the rank of this processor */
    int size = controller->GetNumberOfProcesses();

    // If not master process, do the vtkMarchingCubes implementation
    if (rank != 0)
    {
        const char* prefix = argv[2];
        process(rank, size, controller, prefix);
    }

    // Master
    else
    {   
        // to append each piece into 1 big vtk file
        vtkAppendPolyData *appendWriter = vtkAppendPolyData::New();

        //VT_ON();
        //VT_USER_START("Region 5");

        // go through the processes, and append
        for(int k = 1; k < size; k++)
        {
            vtkPolyData* pd = vtkPolyData::New();
            
            controller->Receive(pd, k, 101);

            #if VTK_MAJOR_VERSION <= 5
                appendWriter->AddInput(pd);
            #else
                appendWriter->AddInputData(pd);
            #endif

            appendWriter->Update();
        }

        vtkPolyDataWriter *pWriter = vtkPolyDataWriter::New();
        pWriter->SetFileName(argv[1]);

        #if VTK_MAJOR_VERSION <= 5
            pWriter->SetInput(appendWriter->GetOutput());
        #else
            pWriter->SetInputData(appendWriter->GetOutputPort());
        #endif

        // output vtk file
        pWriter->Write();

        double t2 = MPI_Wtime();

        double time = t2 - t1;

        printf("MPI_Wtime measured the time elapsed to be: %f\n", time);

        FILE *out_file = fopen("../OutPutFile.txt", "a"); // write only 

        char strPD[6];

        sprintf(strPD, "%f \n", time);

        fputs(strPD, out_file);
        fclose(out_file);
    }

    controller->Finalize(); 
    controller->Delete();

    return EXIT_SUCCESS;
}
