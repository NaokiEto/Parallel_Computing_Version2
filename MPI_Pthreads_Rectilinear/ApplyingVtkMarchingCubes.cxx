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
#include <vtkMPIController.h>

#include "/home/users/neto/NaokiEto/openmpi-install/include/mpi.h"
#include <stdio.h>
#include <pthread.h>
#include "/home/users/neto/NaokiEto/vt-install/include/vampirtrace/vt_user.h"
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

/**
 * This program converts a vtkPolyData image into volume representation (vtkImageData) 
 * where the foreground voxels are 1 and the background voxels are 0. Internally 
 * vtkPolyDataToImageStencil is utilized. The resultant image is saved to disk in 
 * metaimage file formats.
 */

typedef struct Param_Function
{
    int NumThreads;
    const char *VTKinput;
    int ThreadId;
    vtkPolyData* vtkPiece;
    int procRank;
} params;


void* thread_function(void* ptr)
{
    params* NewPtr;
    NewPtr = (params*) ptr;

    /* The vtk file extension we want to search for */

    int NumOfCharPD = strlen(NewPtr->VTKinput) + 5 + 1;
    int original = NumOfCharPD;

    // figure out how many characters are in the temporary file name, 
    if((NewPtr->procRank-1)*(NewPtr->NumThreads) + NewPtr->ThreadId > 9)
        NumOfCharPD += (int) log10((NewPtr->procRank-1)*(NewPtr->NumThreads) + NewPtr->ThreadId);

    char buf[(int) NumOfCharPD - original + 1];
 
    sprintf(buf, "%d", (NewPtr->procRank-1)*(NewPtr->NumThreads) + NewPtr->ThreadId);

    const char* suffix = ".vtk";

    const char* prefix = NewPtr->VTKinput;

    char prefix_suffix[NumOfCharPD];

    strcpy(prefix_suffix, prefix);

    strcat(prefix_suffix, buf);

    strcat(prefix_suffix, suffix);

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

    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    

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

    // Pass the meta image data to the vtkMarchingCubes
    vtkExtractVOI *voi = vtkExtractVOI::New(); 
    voi->SetInput(imgstenc->GetOutput());  
    voi->SetSampleRate(1,1,1);

    voi->Update();

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

    #if VTK_MAJOR_VERSION <= 5
        NewPtr->vtkPiece = triangleCellNormals->GetOutput();
    #else
        NewPtr->vtkPiece = triangleCellNormals->GetOutputPort();  
    #endif
}

int main(int argc, char *argv[])
{
    //VT_ON();

    vtkMPIController* controller = vtkMPIController::New();

    double t1 = MPI_Wtime();

    // Initializing MPI
    controller->Initialize(&argc, &argv);

    /* Figure out total amount of processors */
    int MPI_rank = controller->GetLocalProcessId();

    /* Figure out the rank of this processor */
    int MPI_size = controller->GetNumberOfProcesses();

    /* The master process will be of rank 0 */
    int MASTER = 0;

    int pthreads_size = atoi(argv[1]);

	pthread_t threads[pthreads_size];

    // this is the input into the function for each thread
    params thread_data_array[pthreads_size];

    //VT_OFF();

    // If not master process, do the vtkMarchingCubes implementation
    if (MPI_rank >= 1)
    {
        //creating threads
        // The file extension to look for is .f.vtk
        std::string fileprefix = argv[3];

	    for (int f = 0; f < pthreads_size; f++) {      
            //creating threads
            thread_data_array[f].NumThreads = pthreads_size;
            thread_data_array[f].procRank = MPI_rank;
            thread_data_array[f].VTKinput = fileprefix.c_str();
            thread_data_array[f].ThreadId = f;
		    pthread_create(&threads[f], NULL, thread_function, (void*)&thread_data_array[f]);
        }

	    for (int j = 0; j < pthreads_size; j++)
        {
		    pthread_join(threads[j], NULL);
        }

        vtkAppendPolyData *appendWriter = vtkAppendPolyData::New();

        for(int y = 0; y < pthreads_size; y++)
        {
            vtkPolyData *inputNum = vtkPolyData::New();

            inputNum->ShallowCopy(thread_data_array[y].vtkPiece);

            #if VTK_MAJOR_VERSION <= 5
                appendWriter->AddInput(thread_data_array[y].vtkPiece);
            #else
                appendWriter->AddInputData(inputNum);
            #endif

            appendWriter->Update();
        }

        // send the vtkPolyData to the master process
        #if VTK_MAJOR_VERSION <= 5
            controller->Send(appendWriter->GetOutput(), 0, 1);
        #else
            controller->Send(appendWriter->GetOutputPort(), 0, 1);  
        #endif
    }

    if (MPI_rank == MASTER)
    {
        VT_ON();
        VT_USER_START("FINISH");

        // to append each piece into 1 big vtk file
        vtkAppendPolyData *appendWriterMASTER = vtkAppendPolyData::New();

        // go through the processes, and append
        for(int k = 1; k < MPI_size; k++)
        {
            vtkPolyData* pd = vtkPolyData::New();
            controller->Receive(pd, k, 1);

            #if VTK_MAJOR_VERSION <= 5
                appendWriterMASTER->AddInput(pd);
            #else
                appendWriterMASTER->AddInputData(pd);
            #endif

            appendWriterMASTER->Update();
        }

        vtkPolyDataWriter *pWriter = vtkPolyDataWriter::New();
        pWriter->SetFileName(argv[2]);

        #if VTK_MAJOR_VERSION <= 5
            pWriter->SetInput(appendWriterMASTER->GetOutput());
        #else
            pWriter->SetInputData(appendWriterMASTER->GetOutputPort());
        #endif

        pWriter->Write();

        VT_USER_END("FINISH");
        VT_OFF();

        double t2 = MPI_Wtime();

        double time = t2 - t1;

        printf("MPI_Wtime measured the time elapsed to be: %f\n", time);
/*
        FILE *out_file = fopen("../OutPutFile.txt", "a"); // write only 

        char strPD[6];

        sprintf(strPD, "%f \n", time);

        fprintf(out_file, strPD);
        fclose(out_file);
*/

/*
        // Remove any duplicate points.
        vtkCleanPolyData *cleanFilter = vtkCleanPolyData::New();
        cleanFilter->SetInputConnection(appendWriterMASTER->GetOutputPort());
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
