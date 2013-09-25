/*
 * =====================================================================================
 *
 *       Filename:  Visualization.cxx
 *
 *    Description:  Visualization (End of Pipeline)
 *
 *        Version:  1.0
 *        Created:  09/04/2012 01:48:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

// =======
// Headers
// =======

#include <Visualization.h>
#include <Reader.h>

// For Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// Keys
#include <vtkInformationDoubleVectorKey.h>
#include <vtkInformationKey.h>
#include <FilterInformation.h>

// String
#include <cstring>

// General
#include <vtkSmartPointer.h>

// Data
#include <vtkStructuredGrid.h>
#include <vtkDataSet.h>

// Visualization
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCallbackCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>

// For DEBUG
#include <vtkInformationRequestKey.h>

// ======
// Macros
// ======

vtkStandardNewMacro(Visualization);
vtkCxxRevisionMacro(Visualization,"$Revision 1.0$");

// ===========
// Constructor
// ===========

Visualization::Visualization()
{
    // Member Data
    this->VisualizationInformation = FilterInformation::New();

    // Pipeline
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

// ==========
// Destructor
// ==========

Visualization::~Visualization()
{
    if(this->VisualizationInformation != NULL)
    {
        VisualizationInformation->Delete();
    }
    this->VisualizationInformation = NULL;
}

// ==========
// Print Self
// ==========

void Visualization::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ===================
// Accessors, Mutators
// ===================


// ===============
// Process Request
// ===============

int Visualization::ProcessRequest(
        vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    DEBUG(<< request->GetRequest()->GetName());

    //  Request Data Object
    if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
        return this->RequestDataObject(request,inputVector,outputVector);
    }

    // Request Infomration
    if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
        return this->RequestInformation(request,inputVector,outputVector);
    }

    // Request Update Extent
    if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
        return this->RequestUpdateExtent(request,inputVector,outputVector);
    }

    // Request Data
    if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
        return this->RequestData(request,inputVector,outputVector);
    }

    // Otherwise use superclass
    return this->Superclass::ProcessRequest(request,inputVector,outputVector);
}

// ===========================
// Fill Input Prot Information
// ===========================

int Visualization::FillInputPortInformation(int port,vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkStructuredGrid");
        DEBUG(<< "Success");
        return 1;
    }
    DEBUG(<< "Failure");
    return 0;
}

// ============================
// Fill Output Port Information
// ============================

int Visualization::FillOutputPortInformation(int port, vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkStructuredGrid");
        DEBUG(<< "Success");
        return 1;
    }
    DEBUG(<< "Failure");
    return 0;
}

// ===================
// Request Data Object
// ===================

int Visualization::RequestDataObject(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // Ports
    for(unsigned int i=0; i<static_cast<unsigned int>(this->GetNumberOfOutputPorts()); i++)
    {
        // Output
        vtkInformation *outputInfo = outputVector->GetInformationObject(i);
        vtkStructuredGrid *output = vtkStructuredGrid::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

        // Create new instance
        if(output == NULL)
        {
            output = vtkStructuredGrid::New();
            outputInfo->Set(vtkDataObject::DATA_OBJECT(),output);
            output->FastDelete();
            this->GetOutputPortInformation(i)->Set(vtkDataObject::DATA_EXTENT_TYPE(),output->GetExtentType());
        }
    }

    DEBUG(<< "Success");
    return 1;
}

// ===================
// Request Information
// ===================

int Visualization::RequestInformation(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input 
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // 1- TIME STEPS //

    // Get Time Steps from Input 

    // FilterInformation VisualizationInformation;
    vtkInformationDoubleVectorKey *InputTimeStepsKey = 
        this->VisualizationInformation->GetDoubleVectorKey(inputInfo,FilterInformation::TIME_STEPS());

    double *TimeSteps = inputInfo->Get(InputTimeStepsKey);
    unsigned int TimeStepsLength = inputInfo->Length(InputTimeStepsKey);

    if(TimeStepsLength < 1)
    {
        vtkErrorMacro("TimeStepsLength is zero");
        return 0;
    }

    if(TimeSteps == NULL)
    {
        DEBUG(<< "TIME_STEPS is NULL.")
    }

    DISPLAY(TimeSteps,TimeStepsLength)

    // Set Time Steps to Output
    outputInfo->Set(FilterInformation::TIME_STEPS(),TimeSteps,TimeStepsLength);

    // 2- TIME RANGE //

    // Get Time Range from Input
    vtkInformationDoubleVectorKey *InputTimeRangeKey = 
        this->VisualizationInformation->GetDoubleVectorKey(inputInfo,FilterInformation::TIME_RANGE());
    double *TimeRange = inputInfo->Get(InputTimeRangeKey);

    if(TimeRange == NULL)
    {
        vtkErrorMacro("TimeRange is NULL");
        return 0;
    }

    DISPLAY(TimeRange,2)

    // Set Time Range to Output
    outputInfo->Set(FilterInformation::TIME_RANGE(),TimeRange,2);

    DEBUG(<< "Success");

    return 1;
}

// =====================
// Request Update Extent
// =====================

int Visualization::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    double UpdateTimeSteps[2] = {1,2};
    unsigned int UpdateTimeStepsLength = 2;

    // Set Update Time Steps to Input
    inputInfo->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);

    return 1;
}

// ============
// Request Data
// ============

int Visualization::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkStructuredGrid *input = vtkStructuredGrid::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));

    if(input == NULL)
    {
        ERROR(<< "input is NULL");
        vtkErrorMacro("input is NULL");
    }
    
    if(input->GetNumberOfPoints() < 1)
    {
        ERROR(<< "Number of input points is zero.");
        vtkErrorMacro("Number of input points is zero.");
        return 0;
    }

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkStructuredGrid *output = vtkStructuredGrid::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Visualization //
    // this->VisualizeDataObject();

    return 1;
}

// ====================
// Visualize DataObject
// ====================

void Visualization::VisualizeDataObject(vtkDataObject *DataObject)
{
    // Cast DataObject to StructuredGrid
    vtkSmartPointer<vtkStructuredGrid> StructuredGrid = vtkStructuredGrid::SafeDownCast(DataObject);

    // Mapper
    vtkSmartPointer<vtkDataSetMapper> Mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    Mapper->SetInputData(StructuredGrid);

    // Property
    vtkSmartPointer<vtkProperty> Property = vtkSmartPointer<vtkProperty>::New();

    // Actor
    vtkSmartPointer<vtkActor> Actor = vtkSmartPointer<vtkActor>::New();
    Actor->SetMapper(Mapper);
    Actor->SetProperty(Property);

    // Renderer
    vtkSmartPointer<vtkRenderer> Renderer = vtkSmartPointer<vtkRenderer>::New();
    Renderer->AddActor(Actor);

    // Render Window
    vtkSmartPointer<vtkRenderWindow> RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    RenderWindow->AddRenderer(Renderer);

    // Interactor
    vtkSmartPointer<vtkRenderWindowInteractor> Interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    Interactor->SetRenderWindow(RenderWindow);

    // Start Graphics
    Interactor->Initialize();
    Renderer->Render();
    Interactor->Start();
}
