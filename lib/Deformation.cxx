/*
 * =====================================================================================
 *
 *       Filename:  Deformation.cxx
 *
 *    Description:  Computes Cauchy Green Tensor
 *
 *        Version:  1.0
 *        Created:  11/17/2012 03:36:13 PM
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

#include <Deformation.h>

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// Keys
#include <vtkInformationDoubleVectorKey.h>
#include <FilterInformation.h>

// Progress
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

// Computation
#include <vtkStructuredGrid.h>

// DEBUG
#include <vtkInformationRequestKey.h>

// ======
// Macros
// ======

vtkStandardNewMacro(Deformation);
// vtkCxxRevisionMacro(Deformation,"$Revision 1.0$");

// ===========
// Constructor
// ===========

Deformation::Deformation()
{
    // Member Data


    // Pipeline
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(2);
}

// ==========
// Destructor
// ==========

Deformation::~Deformation()
{
    // TO BE ADDED LATER
}

// ==========
// Print Self
// ==========

void Deformation::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ===================
// Accessors, Mutators
// ===================


// ==========
// Get Output
// ==========

vtkStructuredGrid * Deformation::GetOutput()
{
    return this->GetOutput(0);
}

vtkStructuredGrid * Deformation::GetOutput(int port)
{
    return vtkStructuredGrid::SafeDownCast(this->GetOutputDataObject(port));
}

// ==========
// Set Output
// ==========

void Deformation::SetOutput(vtkDataObject *OutputDataObject)
{
    this->SetOutput(0,OutputDataObject);
}

void Deformation::SetOutput(int port, vtkDataObject *OutputDataObject)
{
    this->GetExecutive()->SetOutputData(port,OutputDataObject);
}

// =========
// Get Input
// =========

vtkDataSet * Deformation::GetInput()
{
    return this->GetInput(0);
}

vtkDataSet * Deformation::GetInput(int port)
{
    return vtkDataSet::SafeDownCast(this->GetExecutive()->GetInputData(port,0));
}

// =========
// Set Input
// =========

void Deformation::SetInput(vtkDataObject *InputDataObject)
{
    this->SetInput(0,InputDataObject);
}

void Deformation::SetInput(int port, vtkDataObject *InputDataObject)
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

void Deformation::AddInput(vtkDataObject *InputDataObject)
{
    this->AddInput(0,InputDataObject);
}

void Deformation::AddInput(int port, vtkDataObject *InputDataObject)
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

int Deformation::ProcessRequest(
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

int Deformation::FillInputPortInformation(int port, vtkInformation *info)
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

int Deformation::FillOutputPortInformation(int port, vtkInformation *info)
{
    switch (port)
    {
        // Output Port 0 is connected to LCS Filter
        case 0:
        {
            info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkStructuredGrid");
            DEBUG(<< "Success");
            return 1;
        }

        // Output Prot 1 is connected to Visualization Filter
        case 1:
        {
            info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkStructuredGrid");
            DEBUG(<< "Success");
            return 1;
        }

        // Output Port is not supported
        default:
        {
            vtkErrorMacro("Port is not suppoerted.");
            DEBUG(<< "Failure");
            return 0;
        }
    }
}

// ===================
// Request Data Object
// ===================

int Deformation::RequestDataObject(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // Ports
    for(unsigned int port=0; port < static_cast<unsigned int>(this->GetNumberOfOutputPorts()); port++)
    {
        // Output Port 0 is connected to LCS Filter
        if(port == 0)
        {
            vtkInformation *outputInfo0 = outputVector->GetInformationObject(port);
            vtkStructuredGrid *output0 = vtkStructuredGrid::SafeDownCast(outputInfo0->Get(vtkDataObject::DATA_OBJECT()));

            // Create new instance
            if(output0 == NULL)
            {
                output0 = vtkStructuredGrid::New();
                outputInfo0->Set(vtkDataObject::DATA_OBJECT(),output0);
                output0->FastDelete();
                this->GetOutputPortInformation(port)->Set(vtkDataObject::DATA_EXTENT_TYPE(),output0->GetExtentType());
            }
        }

        // Output Prot 1 is connected to Visualization Filter
        else if(port == 1)
        {
            vtkInformation *outputInfo1 = outputVector->GetInformationObject(port);
            vtkStructuredGrid *output1 = vtkStructuredGrid::SafeDownCast(outputInfo1->Get(vtkDataObject::DATA_OBJECT()));

            // Create new instance
            if(output1 == NULL)
            {
                output1 = vtkStructuredGrid::New();
                outputInfo1->Set(vtkDataObject::DATA_OBJECT(),output1);
                output1->FastDelete();
                this->GetOutputPortInformation(port)->Set(vtkDataObject::DATA_EXTENT_TYPE(),output1->GetExtentType());
            }
        }
        
        // No Valid output port
        else
        {
            DEBUG(<< "Failure");
            vtkErrorMacro("Port is not supperted.");
            return 0;
        }
    }

    DEBUG(<< "Success");
    return 1;
}

// ===================
// Request Information
// ===================

int Deformation::RequestInformation(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input Info
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Output Info
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // 1- Time Steps //

    // Check if TimeSteps key exists in inputInfo
    if(!inputInfo->Has(FilterInformation::TIME_STEPS()))
    {
        ERROR(<< "inputInfo does not have TIME_STEPS key.");
        vtkErrorMacro("inputInfo does not have TIME_STEPS key.");
        return 0;
    }

    // Get TimeSteps
    unsigned int TimeStepsLength = inputInfo->Length(FilterInformation::TIME_STEPS());
    double * TimeSteps = inputInfo->Get(FilterInformation::TIME_STEPS());

    // Check TimeSteps
    if(TimeStepsLength < 1)
    {
        ERROR(<< "TimeStepsLength is zero.");
        vtkErrorMacro("TimeStepsLength is zero.");
        return 0;
    }

    if(TimeSteps == NULL)
    {
        ERROR(<< "TimeSteps is NULL.");
        vtkErrorMacro("TimeSteps is NULL.");
        return 0;
    }

    // Set TimeSteps to outputInfo
    outputInfo->Set(FilterInformation::TIME_STEPS(),TimeSteps,TimeStepsLength);

    // 2- Time Range //

    // Check if TimeRange Key is in inputInfo
    if(!inputInfo->Has(FilterInformation::TIME_RANGE()))
    {
        ERROR(<< "inputInfo does not have TIME_RANGE.");
        vtkErrorMacro("inputInfo does not have TIME_RANGE.");
        return 0;
    }

    // Get TimeRange
    double *TimeRange = inputInfo->Get(FilterInformation::TIME_RANGE());

    // Check TimeRange
    if(TimeRange == NULL)
    {
        ERROR(<< "TimeRange is NULL.");
        vtkErrorMacro("TimeRange is NULL.");
        return 0;
    }

    // Set TimeRange to outputInfo
    outputInfo->Set(FilterInformation::TIME_RANGE(),TimeRange,2);

    DEBUG(<< "Sucess");
    return 1;
}

// =====================
// Request Update Extent
// =====================

int Deformation::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // input Info
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Check inputInfo
    if(inputInfo == NULL)
    {
        ERROR(<< "inputInfo is NULL.");
        vtkErrorMacro("inputInfo is NULL.");
        return 0;
    }

    // output Info
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // Check outputInfo
    if(outputInfo == NULL)
    {
        ERROR(<< "outputInfo is NULL.");
        vtkErrorMacro("outputInfo is NULL.");
        return 0;
    }

    // Update Time Steps //

    // Check if UpdateTimeSteps Key exists in outputInfo
    if(!outputInfo->Has(FilterInformation::UPDATE_TIME_STEPS()))
    {
        ERROR(<< "outputInfo does not have UPDATE_TIME_STEPS key.");
        vtkErrorMacro("outputInfo does not have UPDATE_TIME_STEPS key.");
        return 0;
    }

    // Get UpdateTimeSteps
    unsigned int UpdateTimeStepsLength = outputInfo->Length(FilterInformation::UPDATE_TIME_STEPS());
    double *UpdateTimeSteps = outputInfo->Get(FilterInformation::UPDATE_TIME_STEPS());

    // Check UpdateTimeSteps
    if(UpdateTimeStepsLength < 1)
    {
        ERROR(<< "UpdatTimeStepsLength is zero.");
        vtkErrorMacro("UpdateTimeStepsLength is zero.");
        return 0;
    }

    if(UpdateTimeSteps == NULL)
    {
        ERROR(<< "UpdateTimeSteps is NULL.");
        vtkErrorMacro("UpdateTimeSteps is NULL.");
        return 0;
    }

    // Set UpdateTimeSteps to inputInfo
    inputInfo->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);

    DEBUG(<< "Success");
    return 1;
}

// ============
// Request Data
// ============

int Deformation::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkStructuredGrid *input = vtkStructuredGrid::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Check input
    if(input == NULL)
    {
        ERROR(<< "input is NULL");
        vtkErrorMacro("input is NULL");
        return 0;
    }
    else if(input->GetNumberOfPoints() == 0)
    {
        ERROR(<< "Number of input points is zero.");
        vtkErrorMacro("Number of input points is zero.");
        return 0;
    }

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkStructuredGrid *output = vtkStructuredGrid::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Computation
    output->ShallowCopy(input);

    DEBUG(<< "Success");
    return 1;
}
