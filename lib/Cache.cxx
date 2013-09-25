/*
 * =====================================================================================
 *
 *       Filename:  Cache.cxx
 *
 *    Description:  Cache
 *
 *        Version:  1.0
 *        Created:  09/11/2012 02:25:50 PM
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

#include <Cache.h>

// For Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// For Keys
#include <vtkInformationDoubleVectorKey.h>
#include <FilterInformation.h>

// For Computation
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>

// For DEBUG
#include <vtkInformationRequestKey.h>

// for Backward Compability
#include <vtkVersion.h>

// ======
// Macros
// ======

vtkStandardNewMacro(Cache);
// vtkCxxRevisionMacro(Cache,"$Revision 1.0$");

// ===========
// Constructor
// ===========

Cache::Cache()
{
    // Member Data
    this->CacheSize = 2;

    // Pipeline
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

// ==========
// Destructor
// ==========

Cache::~Cache()
{
}

// ==========
// Print Self
// ==========

void Cache::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ==========
// Get Output
// ==========

vtkMultiBlockDataSet * Cache::GetOutput()
{
    return this->GetOutput(0);
}

vtkMultiBlockDataSet * Cache::GetOutput(int port)
{
    return vtkMultiBlockDataSet::SafeDownCast(this->GetOutputDataObject(port));
}

// ==========
// Set Output
// ==========

void Cache::SetOutput(vtkDataObject *OutputDataObject)
{
    this->SetOutput(0,OutputDataObject);
}

void Cache::SetOutput(int port, vtkDataObject *OutputDataObject)
{
    this->GetExecutive()->SetOutputData(port,OutputDataObject);
}

// =========
// Get Input
// =========

vtkDataSet * Cache::GetInput()
{
    return this->GetInput(0);
}

vtkDataSet * Cache::GetInput(int port)
{
    return vtkDataSet::SafeDownCast(this->GetExecutive()->GetInputData(port,0));
}

// =========
// Set Input
// =========

void Cache::SetInput(vtkDataObject *InputDataObject)
{
    this->SetInput(0,InputDataObject);
}

void Cache::SetInput(int port, vtkDataObject *InputDataObject)
{
    if(InputDataObject != NULL)
    {
        #if VTK_MAJOR_VERSION <= 5
        this->SetInputConnection(port,InputDataObject->GetProducerPort());
        #else
        this->SetInputData(port,InputDataObject);
        #endif
    }
    else
    {
        this->SetInputConnection(port,0);
    }
}

// =========
// Add Input
// =========

void Cache::AddInput(vtkDataObject *InputDataObject)
{
    this->AddInput(0,InputDataObject);
}

void Cache::AddInput(int port, vtkDataObject *InputDataObject)
{
    if(InputDataObject != NULL)
    {
        #if VTK_MAJOR_VERSION <= 5
        this->AddInputConnection(port,InputDataObject->GetProducerPort());
        #else
        this->AddInputDataObject(port,InputDataObject);
        #endif
    }
    else
    {
        this->AddInputConnection(port,0);
    }
}

// ===============
// Process Request
// ===============

int Cache::ProcessRequest(
        vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    DEBUG(<< request->GetRequest()->GetName());

    // Request Data Object
    if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
        return this->RequestDataObject(request,inputVector,outputVector);
    }

    // Request Information
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
// Fill Input Port Information
// ===========================

int Cache::FillInputPortInformation(int port,vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkMultiBlockDataSet");
        DEBUG(<< "Success");
        return 1;
    }

    DEBUG(<< "Failure");
    return 0;
}

// ============================
// Fill Output Port Information
// ============================

int Cache::FillOutputPortInformation(int port, vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkMultiBlockDataSet");
        DEBUG(<< "Success");
        return 1;
    }

    DEBUG(<< "Failure");
    return 0;
}

// ===================
// Request Data Object
// ===================

int Cache::RequestDataObject(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // Ports
    for(unsigned int i=0; i<static_cast<unsigned int>(this->GetNumberOfOutputPorts()); i++)
    {
        // Output
        vtkInformation *outputInfo = outputVector->GetInformationObject(i);
        vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

        // Create new instance
        if(output == NULL)
        {
            output = vtkMultiBlockDataSet::New();
            outputInfo->Set(vtkDataObject::DATA_OBJECT(),output);
            output->FastDelete();
            #if VTK_MAJOR_VERSION <= 5
            output->SetPipelineInformation(outputInfo);
            #endif
            this->GetOutputPortInformation(i)->Set(vtkDataObject::DATA_EXTENT_TYPE(),output->GetExtentType());
        }
    }

    DEBUG(<< "Success");
    return 1;
}

// ===================
// Request Information
// ===================

int Cache::RequestInformation(
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
    unsigned int TimeStepsLength = inputInfo->Length(FilterInformation::TIME_STEPS());
    double *TimeSteps = inputInfo->Get(FilterInformation::TIME_STEPS());

    // Set Time Steps to Output
    outputInfo->Set(FilterInformation::TIME_STEPS(),TimeSteps,TimeStepsLength);

    // 2- TIME RANGE //

    // Get Time Range from Input
    double *TimeRange = inputInfo->Get(FilterInformation::TIME_RANGE());

    // Set Time Range to Output
    outputInfo->Set(FilterInformation::TIME_RANGE(),TimeRange,2);

    // 3- UPDATE TIME STEPS //

    // Get Update Time Steps from Input
    unsigned int UpdateTimeStepsLength = inputInfo->Length(FilterInformation::UPDATE_TIME_STEPS());
    double *UpdateTimeSteps = inputInfo->Get(FilterInformation::UPDATE_TIME_STEPS());

    // Set Update Time Steps to Output
    outputInfo->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);

    DISPLAY(TimeSteps,TimeStepsLength);
    DISPLAY(TimeRange,2);
    DISPLAY(UpdateTimeSteps,UpdateTimeStepsLength);
    DEBUG(<< "Success");
    return 1;
}

// =====================
// Request Update Extent
// =====================

int Cache::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    DEBUG(<< "Pass");
    return 1;
}

// ============
// Request Data
// ============

int Cache::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet *input = vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Computations
    output->ShallowCopy(input);

    DEBUG(<< "Pass");
    return 1;
}
