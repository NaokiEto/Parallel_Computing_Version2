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

//#include "/work/naoki/openmpi-1.6.5/ompi/contrib/vt/vt/vtlib/vt_pthreadreg.h"

__VT_EXTERN_DECL int VT_pthread_create__(pthread_t* thread, const pthread_attr_t* attr,
                               void *(*start_routine)(void*), void* arg);
__VT_EXTERN_DECL int VT_pthread_join__(pthread_t thread, void** value_ptr);

#define pthread_create VT_pthread_create__
#define pthread_join VT_pthread_join__

// holds search results
std::vector<std::string> results;

/**
* This recursive search algorithm function will be used later to find the vtk 
* file in the directory. The user is to place the vtk file in the build 
* directory, and this function will find it and output it.
*/
void search(std::string curr_directory, std::string extension){
	DIR* dir_point = opendir(curr_directory.c_str());
	dirent* entry = readdir(dir_point);

	while (entry){
	    // filename
		std::string fname = entry->d_name;
		// if filename's last characters are extension
		if (fname.find(extension, (fname.length() - extension.length())) != std::string::npos)
            // add filename to results vector
			results.push_back(fname);

		entry = readdir(dir_point);
	}
	return;
}

/**
 * This struct contains the id of the thread, the filename of the vtk file 
 * that was found, and the number of total threads. This will be useful
 * for determining the bounds in the vtk data for each particular thread.
*/
typedef struct Param_Function
{
    const char *VTKInput;
    int threadId;
    int numThreads;
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
    VT_USER_START("Region 2");

    // The vtkPolyDataReader to read the input vtk file in the build
    // directory
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(NewPtr->VTKInput);
    reader->Update();

    VT_USER_END("Region 2");
    VT_OFF();

    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    

    reader->Update();

    VT_ON();
    VT_USER_START("Region 3"); 

    // desired volume spacing
    double spacing[3];
    spacing[0] = 0.1;
    spacing[1] = 0.1;
    spacing[2] = 0.1;
    whiteImage->SetSpacing(spacing);

    /* The bounds that this processor will deal with in the vtk file */
    double bounds[6];
    vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();

    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
    reader->Update();

    // Each pthread will work with a different piece of vtk data. It's basically just 
    // a rectangular prism of data
    bounds[0] = -10.0 + 20.0 * double(NewPtr->threadId) / double(NewPtr->numThreads);
    bounds[1] = -10.0 + 20.0 * double(NewPtr->threadId + 1.0) / double(NewPtr->numThreads);
    bounds[2] = -10;
    bounds[3] = 10;
    bounds[4] = -10;
    bounds[5] = 10;

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
    origin[0] = bounds[0] + spacing[0]/10;
    origin[1] = bounds[2] + spacing[1]/10;
    origin[2] = bounds[4] + spacing[2]/10;
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

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();

    #if VTK_MAJOR_VERSION <= 5
        pol2stenc->SetInput(reader->GetOutput());
    #else
        pol2stenc->SetInputData(reader->GetOutput());
    #endif

    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
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
    int index = 31;

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
    int NumOfCharPD = 18;

    // figure out how many characters are in the temporary file name, 
    // i.e. "ShrimpChowFun128475.vtk"
    if(NewPtr->threadId >= 1)
        NumOfCharPD += log10(NewPtr->threadId);

    char strPD[NumOfCharPD];

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
    /* Number of threads */
    int size = atoi(argv[1]);

    /* Array with elements of type params (the structure defined above) */
    params thread_data_array[size];

    /* for the joining of threads later on */
    void* exit_status;

    // The file extension to look for is vtk
    std::string extension;
    extension = "vtk";

    VT_ON();
    VT_USER_START("Region 1");

	// setup search parameters
	std::string curr_directory = get_current_dir_name();

	search(curr_directory, extension);
    std::string result = results[0];

    VT_USER_END("Region 1");
    VT_OFF();

	pthread_t threads[size];

	for (int i = 0; i < size; i++) {      
        // Creating threads
	    thread_data_array[i].VTKInput = result.c_str();
        thread_data_array[i].numThreads = size;
        thread_data_array[i].threadId = i;
		pthread_create(&threads[i], NULL, thread_function, (void*)&thread_data_array[i]);
	}

    VT_ON();
    VT_USER_START("Region 5");

    // Joining threads to make sure all threads are terminated
	for (int j = 0; j < size; j++)
    {
		pthread_join(threads[j], NULL);
    }

    /* to append each piece into 1 big vtk file */
    vtkAppendPolyData *appendWriter = vtkAppendPolyData::New();

    // Go through the temporary files and append the vtk data
    for(int k = 0; k < size; k++)
    {
        vtkPolyData *inputNum = vtkPolyData::New();
        vtkPolyDataReader *reader = vtkPolyDataReader::New();

        int NumOfCharPDMASTER;

        NumOfCharPDMASTER = 18;

        /* number of characters in file name */
        if (k > 0)
            NumOfCharPDMASTER = log10(k) + 18;

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
	return EXIT_SUCCESS;
}
