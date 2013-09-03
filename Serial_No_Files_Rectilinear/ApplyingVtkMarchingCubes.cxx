/**
* Do whatever you want with public license
* Version 1, August 13, 2013
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
* 1. Uses VTK_MAJOR_VERSION <= 5
* 
*/
/**
* @file ApplyingVtkMarchingCubes.cxx
* @author Naoki Eto
* @date August 13, 2013
* @brief This program gets the VTK file and converts it to a metaimage data 
*        so that vtkMarchingCubes class can be applied. It applies marching 
*        cubes to the data and outputs the resulting vtk file.
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
#include <vtkPointData.h>
#include <string.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataArray.h>
#include <vtkLookupTable.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>

#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridReader.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkContourFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkPoints.h>
#include <algorithm> 

#include <time.h>

#include <stdio.h>
#include "/work2/vt-system-install/include/vampirtrace/vt_user.h"

/**
 * This program converts a vtkPolyData image into volume representation 
 * (vtkImageData) where the foreground voxels are 1 and the background 
 * voxels are 0. Internally vtkPolyDataToImageStencil is utilized. The 
 * resultant image is saved to disk in a metaimage file format. 
 * vtkMarchingCubes is applied to this file format, and a vtk file is 
 * outputted.
 */
int main(int argc, char *argv[])
{
    /* The vtk file extension we want to search for */

    //VT_ON();
    //VT_USER_START("Region 1");

    const char* fp = argv[2];

    int NumOfCharPD = strlen(fp) + 5 + 1;
    int original = NumOfCharPD;

    char buf[(int) NumOfCharPD - original + 1];
 
    sprintf(buf, "%d", 2-1);

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
    reader->Delete();

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

    contour->Update();

    // calc cell normal
    vtkPolyDataNormals *triangleCellNormals= vtkPolyDataNormals::New();

    triangleCellNormals->SetInputConnection(contour->GetOutputPort());

    triangleCellNormals->ComputeCellNormalsOn();
    triangleCellNormals->ComputePointNormalsOff();
    triangleCellNormals->ConsistencyOn();
    triangleCellNormals->AutoOrientNormalsOn();
    triangleCellNormals->Update(); // creates vtkPolyData

    VT_ON();
    VT_USER_START("Region 5");

    /* vtkPolyDataWriter for output vtk file */
    vtkPolyDataWriter *pWriter = vtkPolyDataWriter::New();
    pWriter->SetFileName(argv[1]);

    pWriter->SetInput(triangleCellNormals->GetOutput());

    pWriter->Write();

    triangleCellNormals->Delete();
    pointdata->Delete();
    contour->Delete();

    VT_USER_END("Region 5");
    VT_OFF();
/*
    // Remove any duplicate points.
    vtkCleanPolyData *cleanFilter = vtkCleanPolyData::New();
    cleanFilter->SetInputConnection(triangleCellNormals->GetOutputPort());
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

    //VT_USER_END("ApplyingVtkMarchingCubes");  
*/
    return EXIT_SUCCESS;
}
