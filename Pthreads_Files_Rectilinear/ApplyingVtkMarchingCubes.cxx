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
* @brief This program gets the VTK file,  divides up it up, and performs 
*        pthreading (or POSIX threading). It converts each piece to metaimage 
*        data so that vtkMarchingCubes class can be applied, applies marching 
*        cubes to each process, outputs the resulting vtk data as temporary 
*        files, and then conglomerates the data in the temporary files into 
*        1 vtk file.
* @param[in] argv[1] - number of threads for pthreading (look at 
*            README for more information)
* @param[in] argv[2] - the output's filename
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
#include <vtkMetaImageWriter.h>
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
#include <stdio.h>
#include <pthread.h>
#include <time.h>
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

//#include "/work/naoki/openmpi-1.6.5/ompi/contrib/vt/vt/vtlib/vt_pthreadreg.h"

/**
 * This struct contains the id of the thread and the unique filename of the 
 * vtk file that was found. This will be useful for determining the bounds 
 * in the vtk data for each particular thread.
*/
typedef struct Param_Function
{
    const char * VTKinput;
    int threadId;
} params;

/**
 * This program converts a vtkPolyData image into volume representation 
 * (vtkImageData) where the foreground voxels are 1 and the background 
 * voxels are 0. Internally vtkPolyDataToImageStencil is utilized as 
 * well as pthreads. The resultant image is saved to disk in metaimage 
 * file formats. vtkMarchingCubes is applied to these file formats, and 
 * temporary vtk files are outputted. 
*/
void* thread_function(void* ptr)
{
    /* This is a struct variable that will be useful later on for determining
       the bounds in the vtkpolydata for each thread. */
    params* NewPtr;
    NewPtr = (params*) ptr;

    VT_ON();
    VT_USER_START("Region 1");

    int NumOfCharPD = strlen(NewPtr->VTKinput) + 5 + 1;
    int original = NumOfCharPD;

    // figure out how many characters are in the temporary file name, 
    if(NewPtr->threadId > 9)
        NumOfCharPD += (int) log10(NewPtr->threadId);

    char buf[(int) NumOfCharPD - original + 1];

    sprintf(buf, "%d", NewPtr->threadId);

    const char* suffix = ".vtk";

    const char* prefix = NewPtr->VTKinput;

    char prefix_suffix[NumOfCharPD];

    strcpy(prefix_suffix, prefix);

    strcat(prefix_suffix, buf);

    strcat(prefix_suffix, suffix);

    VT_USER_END("Region 1");
    VT_OFF();

    vtkRectilinearGridReader *reader = vtkRectilinearGridReader::New();

    VT_ON();
    VT_USER_START("Region 2");
    
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

    vtkPolyData* polydata = surfaceFilter->GetOutput();

    /* compute bounds */
    double bounds[6];

    VT_ON();
    VT_USER_START("Region 3");   

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
    for (vtkIdType s = 0; s < count; ++s)
    {
        whiteImage->GetPointData()->GetScalars()->SetTuple1(s, inval);
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

    /* vtkPolyDataWriter for temporary file */
    vtkPolyDataWriter *PDwriter = vtkPolyDataWriter::New();

    /* This is for the temporary file names created. It contains the number 
       of characters in the temporary filename "ShrimpChowFun#.vtk" */
    int NumOfCharPD2 = 18;

    // figure out how many characters are in the temporary file name, 
    // i.e. "ShrimpChowFun128475.vtk"
    if(NewPtr->threadId > 0)
        NumOfCharPD2 = 18 + (int) log10(NewPtr->threadId);

    char strPD[NumOfCharPD2];

    sprintf(strPD, "ShrimpChowFun%d.vtk", NewPtr->threadId);

    PDwriter->SetFileName(strPD);

    #if VTK_MAJOR_VERSION <= 5
	    PDwriter->SetInput(triangleCellNormals->GetOutput());
    #else
	    PDwriter->SetInputData(triangleCellNormals->GetOutputPort());  
    #endif

    PDwriter->Write();

    VT_USER_END("Region 4");
    VT_OFF();
}

/**
 * The main function calls the thread function. The vtkpolydata outputted 
 * from these functions are collected and appended together to form one
 * output vtkpolydata file.
 */
int main(int argc, char *argv[])
{
    struct timespec t0,t1;

    clock_gettime(CLOCK_REALTIME,&t0);

    /* Number of threads */
    int size = atoi(argv[1]);

    /* Array with elements of type params (the structure defined above) */
    params thread_data_array[size];

	pthread_t threads[size];

	for (int f = 0; f < size; f++) {      
        // Creating threads
        // The file extension to look for is .f.vtk
        std::string fileprefix = argv[3];

        thread_data_array[f].VTKinput = fileprefix.c_str();
        thread_data_array[f].threadId = f;
        //pthread_mutex_lock(&lock);
		pthread_create(&threads[f], NULL, thread_function, (void*)&thread_data_array[f]);
        //pthread_mutex_unlock(&lock);
	}

    // Joining threads to make sure all threads are terminated
	for (int j = 0; j < size; j++)
    {
		pthread_join(threads[j], NULL);
    }

    VT_ON();
    VT_USER_START("Region 5");

    /* to append each piece into 1 big vtk file */
    vtkAppendPolyData *appendWriter = vtkAppendPolyData::New();

    // Go through the temporary files and append the vtk data
    for(int k = 0; k < size; k++)
    {
        vtkPolyData *inputNum = vtkPolyData::New();
        vtkPolyDataReader *readerPD = vtkPolyDataReader::New();

        int NumOfCharPDMASTER = 18;

        /* number of characters in file name */
        if (k > 0)
            NumOfCharPDMASTER = (int) log10(k) + 18;

        char strPDMASTER[NumOfCharPDMASTER];

        sprintf(strPDMASTER, "ShrimpChowFun%d.vtk", k);
        readerPD->SetFileName(strPDMASTER);
        readerPD->Update();
        inputNum->ShallowCopy(readerPD->GetOutput());

        #if VTK_MAJOR_VERSION <= 5
            appendWriter->AddInput(readerPD->GetOutput());
        #else
            appendWriter->AddInputData(inputNum);
        #endif

        // remove temporary files
        remove(strPDMASTER);

        appendWriter->Update();
    }

    vtkPolyDataWriter *pWriter = vtkPolyDataWriter::New();
    
    // Output vtkpolydata file
    pWriter->SetFileName(argv[2]);

    #if VTK_MAJOR_VERSION <= 5
        pWriter->SetInput(appendWriter->GetOutput());
    #else
        pWriter->SetInputData(appendWriter->GetOutputPort());
    #endif

    pWriter->Write();

    VT_USER_END("Region 5");
    VT_OFF();

    clock_gettime(CLOCK_REALTIME,&t1);

    double dt = (double) (t1.tv_sec - t0.tv_sec) + (double) (t1.tv_nsec - t0.tv_nsec) / 1.0e9;

    printf("The measured time elapsed to be: %f\n", dt);
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
	//return EXIT_SUCCESS;
    exit(0);
}
