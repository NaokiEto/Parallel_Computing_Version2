/**
* Do whatever you want with public license
* Version 1, August 12, 2013
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
*        the resulting vtk data as temporary files, and then the master 
*        processor (with rank 0) conglomerates the data in the temporary files 
*        into 1 vtk file.
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

#include "/work2/vt-system-install/include/vampirtrace/vt_user.h"
#include <mpi.h>
#include <stdio.h>
#include <math.h>
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

#include <time.h>

/**
 * This program converts a vtkPolyData image into volume representation 
 * (vtkImageData) where the foreground voxels are 1 and the background 
 * voxels are 0. Internally vtkPolyDataToImageStencil is utilized as well
 * as MPI. The resultant image is saved to disk in metaimage file formats. 
 * vtkMarchingCubes is applied to these file formats, and temporary vtk
 * files are outputted. These temporary files are then conglomerated by 
 * the master process, and an output vtk file is outputted.
 */
int main(int argc, char *argv[])
{
    /* Number of processes */
    static int size;

    /* Rank of process */
    int rank;

    MPI_Status status;

    // Initializing MPI
    MPI_Init (&argc, &argv);

    double t1 = MPI_Wtime();

    // Figure out total amount of processors
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Figure out the rank of this processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int NumOfCharPD;

    char strPD[NumOfCharPD];

    // If not master process, do the vtkMarchingCubes implementation
    if (rank >= 1)
    {
        /* This is for the temporary file names created. It contains the number 
           of characters in the temporary filename "ShrimpChowFun#.vtk" */
        NumOfCharPD = 18 + log10(rank);

        /* The vtk file extension we want to search for */

        VT_ON();
        VT_USER_START("Region 1");

        int NumOfCharIN = strlen(argv[2]) + 5 + 1;
        int original = NumOfCharIN;

        // figure out how many characters are in the temporary file name, 
        if(rank - 1 > 9)
            NumOfCharIN += (int) log10(rank - 1);

        char buf[(int) NumOfCharIN - original + 1];
     
        sprintf(buf, "%d", rank-1);

        const char* suffix = ".vtk";

        const char* prefix = argv[2];

        char prefix_suffix[NumOfCharIN];

        strcpy(prefix_suffix, prefix);

        strcat(prefix_suffix, buf);

        strcat(prefix_suffix, suffix);

        VT_USER_END("Region 1");
        VT_OFF();

        VT_ON();
        VT_USER_START("Region 2");

        vtkRectilinearGridReader *reader = vtkRectilinearGridReader::New();

        reader->SetFileName(prefix_suffix);
        reader->Update();


        // Create a grid
        vtkSmartPointer<vtkRectilinearGrid> grid = reader->GetOutput();

        vtkSmartPointer<vtkRectilinearGridToTetrahedra> rectilinearGridToTetrahedra =
        vtkSmartPointer<vtkRectilinearGridToTetrahedra>::New();
        #if VTK_MAJOR_VERSION <= 5
            rectilinearGridToTetrahedra->SetInputConnection(grid->GetProducerPort());
        #else
            rectilinearGridToTetrahedra->SetInputData(grid);
        #endif
            rectilinearGridToTetrahedra->Update();

        vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();

        // Poly Data
        #if VTK_MAJOR_VERSION <= 5
            surfaceFilter->SetInputConnection(rectilinearGridToTetrahedra->GetOutputPort());
        #else
            surfaceFilter->SetInputData(rectilinearGridToTetrahedra->GetOutput());
        #endif
        surfaceFilter->Update(); 

        VT_USER_END("Region 2");
        VT_OFF();

        vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    

       VT_ON();
        VT_USER_START("Region 3");     

        vtkPolyData* polydata = surfaceFilter->GetOutput();

         /* compute bounds */ 
        double bounds[6];

        polydata->GetBounds(bounds);

        /* desired volume spacing */
        double spacing[3];
        spacing[0] = 0.1;
        spacing[1] = 0.1;
        spacing[2] = 0.1;
        whiteImage->SetSpacing(spacing);

        // cut the corresponding white image and set the background:
        vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();

        // polygonal data --> image stencil:
        vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();

        #if VTK_MAJOR_VERSION <= 5
            pol2stenc->SetInput(surfaceFilter->GetOutput());
        #else
            pol2stenc->SetInputData(surfaceFilter->GetOutput());
        #endif

        pol2stenc->SetOutputSpacing(spacing);

        /* compute dimensions */
        int dim[3];
        for (int i = 0; i < 3; i++)
        {
            dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
        }
        whiteImage->SetDimensions(dim);
        whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

        /* compute the origin */
        double origin[3];
        for (int t = 0; t < 3; t++)
        {
            if (bounds[2*t] < 0)
                origin[t] = bounds[2*t] - spacing[t]/10;
            else
                origin[t] = bounds[2*t] + spacing[t]/10;
        }

        whiteImage->SetOrigin(origin);

        #if VTK_MAJOR_VERSION <= 5
	        whiteImage->SetScalarTypeToUnsignedChar();
	        whiteImage->AllocateScalars();
        #else
	       whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
        #endif

        // fill the image with foreground voxels:
        unsigned char inval = 255;
        unsigned char outval = 0;
        vtkIdType count = whiteImage->GetNumberOfPoints();
        for (vtkIdType i = 0; i < count; ++i)
        {
	        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
        }

        pol2stenc->SetOutputOrigin(origin);
        pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
        pol2stenc->Update();

        #if VTK_MAJOR_VERSION <= 5
	        imgstenc->SetInput(whiteImage);
	        imgstenc->SetStencil(pol2stenc->GetOutput());
        #else
	        imgstenc->SetInputData(whiteImage);
	        imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
        #endif

        imgstenc->ReverseStencilOff();
        imgstenc->SetBackgroundValue(outval);
        imgstenc->Update();

        VT_USER_END("Region 3");
        VT_OFF();

        // Pass the meta image data to the vtkMarchingCubes
        vtkExtractVOI *voi = vtkExtractVOI::New(); 
        voi->SetInput(imgstenc->GetOutput());  
        voi->SetSampleRate(1,1,1);

        voi->Update();

        VT_ON();
        VT_USER_START("Region 4");

        //Prepare surface generation
        vtkMarchingCubes *contour = vtkMarchingCubes::New(); //for label images

        contour->SetInput( voi->GetOutput() );

        // choose one label
        int index = 131;

        contour->SetValue(0, index);
        contour->Update(); //needed for GetNumberOfPolys() !!!

        vtkWindowedSincPolyDataFilter *smoother= vtkWindowedSincPolyDataFilter::New();

        #if VTK_MAJOR_VERSION <= 5
	        smoother->SetInput(contour->GetOutput());
        #else
	        smoother->SetInputConnection(contour->GetOutputPort());   
        #endif

        smoother->SetNumberOfIterations(30); // this has little effect on the error!
        smoother->NonManifoldSmoothingOn();
        smoother->NormalizeCoordinatesOn();
        smoother->GenerateErrorScalarsOn();
        smoother->Update();

        vtkPolyData *smoothed_polys = smoother->GetOutput();

        // calc cell normal
        vtkPolyDataNormals *triangleCellNormals= vtkPolyDataNormals::New();

        #if VTK_MAJOR_VERSION <= 5
	        triangleCellNormals->SetInput(smoothed_polys);
        #else
	        triangleCellNormals->SetInputData(smoothed_polys);   
        #endif

        triangleCellNormals->ComputeCellNormalsOn();
        triangleCellNormals->ComputePointNormalsOff();
        triangleCellNormals->ConsistencyOn();
        triangleCellNormals->AutoOrientNormalsOn();
        triangleCellNormals->Update(); // creates vtkPolyData

        vtkPolyDataMapper *mapper= vtkPolyDataMapper::New();

        /* vtkPolyDataWriter for temporary file */
        vtkPolyDataWriter *PDwriter = vtkPolyDataWriter::New();

        sprintf(strPD, "ShrimpChowFun%d.vtk", rank);

        PDwriter->SetFileName(strPD);

        #if VTK_MAJOR_VERSION <= 5
    	    PDwriter->SetInput(triangleCellNormals->GetOutput());
        #else
    	    PDwriter->SetInputData(triangleCellNormals->GetOutputPort());  
        #endif

        PDwriter->Write();

        // Send this temporary file name to the master process
        MPI_Send(strPD, NumOfCharPD, MPI_CHAR, 0, 1, MPI_COMM_WORLD);

        VT_USER_END("Region 4");
        VT_OFF();
    }

    // Master
    else
    {
        /* to append each piece into 1 big vtk file */
        vtkAppendPolyData *appendWriter = vtkAppendPolyData::New();

        VT_ON();
        VT_USER_START("Region 5");

        // Receive all the temporary file names
        for(int i = 1; i < size; i++)
        {
            NumOfCharPD = log10(i) + 18;

            MPI_Recv(strPD, NumOfCharPD, MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
        }

        // Go through the temporary files and append the vtk data
        for(int k = 1; k < size; k++)
        {
            vtkPolyData *inputNum = vtkPolyData::New();
            vtkPolyDataReader *reader = vtkPolyDataReader::New();

            /* number of characters in file name */
            int NumOfCharPDMASTER = log10(k) + 18;

            char strPDMASTER[NumOfCharPDMASTER];

            sprintf(strPDMASTER, "ShrimpChowFun%d.vtk", k);
            reader->SetFileName(strPDMASTER);
            reader->Update();
            inputNum->ShallowCopy(reader->GetOutput());

            #if VTK_MAJOR_VERSION <= 5
                appendWriter->AddInput(reader->GetOutput());
            #else
                appendWriter->AddInputData(inputNum);
            #endif

            // remove temporary files
            remove(strPDMASTER);

            appendWriter->Update();
        }

        // output vtkpolydata file
        vtkPolyDataWriter *pWriter = vtkPolyDataWriter::New();
        pWriter->SetFileName(argv[1]);

        #if VTK_MAJOR_VERSION <= 5
            pWriter->SetInput(appendWriter->GetOutput());
        #else
            pWriter->SetInputData(appendWriter->GetOutputPort());
        #endif

        pWriter->Write();

        VT_USER_END("Region 5");
        VT_OFF();

        double t2 = MPI_Wtime();

        double time = t2 - t1;

        printf("MPI_Wtime measured the time elapsed to be: %f\n", time);

        FILE *out_file = fopen("../OutPutFile.txt", "a"); // write only 

        char strPD[6];

        sprintf(strPD, "%f \n", time);

        fprintf(out_file, strPD);
        fclose(out_file);
/*
        // Remove any duplicate points.
        vtkCleanPolyData *cleanFilter = vtkCleanPolyData::New();
        cleanFilter->SetInputConnection(appendWriter->GetOutputPort());
        cleanFilter->Update();

        // Create a mapper
        vtkPolyDataMapper *cumulativeMapper = vtkPolyDataMapper::New();
        cumulativeMapper->SetInputConnection(cleanFilter->GetOutputPort());

        // Create an actor
        vtkActor *actor = vtkActor::New();
        actor->SetMapper(cumulativeMapper);

        // Create a renderer, render window, and interactor
        vtkRenderer* renderer = vtkRenderer::New();
        vtkRenderWindow *renderWindow = vtkRenderWindow::New();
        vtkRenderWindowInteractor *renderWindowInteractor = vtkRenderWindowInteractor::New();

        // Add the actors to the scene
        renderWindow->AddRenderer(renderer);
        renderWindowInteractor->SetRenderWindow(renderWindow);
        renderer->AddActor(actor);

        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
*/
    }

    MPI_Finalize();    

    return EXIT_SUCCESS;
}
