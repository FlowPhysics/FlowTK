/*
 * =====================================================================================
 *
 *       Filename:  LCS.cxx
 *
 *    Description:  Lagrangian Coherent Structures
 *
 *        Version:  1.0
 *        Created:  11/20/2012 06:31:33 PM
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

#include <LCS.h>

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// Keys
#include <vtkInformationDoubleVectorKey.h>
#include <FilterInformation.h>

// Data
#include <vtkStructuredGrid.h>

// DEBUG
#include <vtkInformationRequestKey.h>


// ======
// Macros
// ======

vtkStandardNewMacro(LCS);
vtkCxxRevisionMacro(LCS,"4Revision 1.0$");

// ===========
// Constructor
// ===========

LCS::LCS()
{
    // Pipeline
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

// ===========
//  Destructor
// ===========

LCS::~LCS()
{
    // TO BE ADDED LATER
}

// ==========
// Print Self
// ==========

void LCS::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ==========
// Get Output
// ==========

vtkStructuredGrid * LCS::GetOutput()
{
    return this->GetOutput(0);
}

vtkStructuredGrid * LCS::GetOutput(int port)
{
    return vtkStructuredGrid::SafeDownCast(this->GetOutputDataObject(port));
}

// ==========
// Set Output
// ==========

void LCS::SetOutput(vtkDataObject *OutputDataObject)
{
    this->SetOutput(0,OutputDataObject);
}

void LCS::SetOutput(int port, vtkDataObject *OutputDataObject)
{
    this->GetExecutive()->SetOutputData(port, OutputDataObject);
}

// =========
// Get Input
// =========

vtkDataSet * LCS::GetInput()
{
    return this->GetInput(0);
}

vtkDataSet * LCS::GetInput(int port)
{
    return vtkDataSet::SafeDownCast(this->GetExecutive()->GetInputData(port,0));
}

// =========
// Set Input
// =========

void LCS::SetInput(vtkDataObject *InputDataObject)
{
    this->SetInput(0,InputDataObject);
}

void LCS::SetInput(int port, vtkDataObject *InputDataObject)
{
    if(InputDataObject != NULL)
    {
        #if VTK_MAJOR_VERSION <= 5
        this->SetInputConnection(port, InputDataObject->GetProducerPort());
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

void LCS::AddInput(vtkDataObject *InputDataObject)
{
    this->AddInput(0,InputDataObject);
}

void LCS::AddInput(int port, vtkDataObject *InputDataObject)
{
    if(InputDataObject != NULL)
    {
        #if VTK_MAJOR_VERSION <= 5
        this->AddInputConnection(port, InputDataObject->GetProducerPort());
        #else
        this->AddInputDataObject(port,InputDataObject);
        #endif
    }
    else
    {
        this->AddInputConnection(port,0);
    }
}

// ===================
// Accessors, Mutators
// ===================


// ===============
// Process Request
// ===============

int LCS::ProcessRequest(
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

    // Otherwise upse superclass
    return this->Superclass::ProcessRequest(request,inputVector,outputVector);
}

// ===========================
// Fill Input Port Information
// ===========================

int LCS::FillInputPortInformation(int port, vtkInformation *info)
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

int LCS::FillOutputPortInformation(int port, vtkInformation *info)
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

int LCS::RequestDataObject(
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

       // Creat new instance
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

int LCS::RequestInformation(
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
    double *TimeSteps = inputInfo->Get(FilterInformation::TIME_STEPS());

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

    // Check is TimeRange key exists in inputInfo
    if(!inputInfo->Has(FilterInformation::TIME_RANGE()))
    {
        ERROR(<< "inputInfo does not have TIME_RANGE key.");
        vtkErrorMacro("inputInfo does not have TIME_RANGE key.");
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

    return 1;
}

// =====================
// Request Update Extent
// =====================

int LCS::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Set Extact Extent //
    unsigned int NumberOfInputPorts = this->GetNumberOfInputPorts();
    for(unsigned int i=0; i<NumberOfInputPorts; i++)
    {
        unsigned int NumberOfInputConnections = this->GetNumberOfInputConnections(i);
        for(unsigned int j=0; j<NumberOfInputConnections; j++)
        {
            vtkInformation *inputInfo = inputVector[i]->GetInformationObject(j);
            inputInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(),1);
        }
    }

    // Input Info
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Output Info
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // Update Time Steps //

    // Check if UpdateTimeSteps key exists in outputInfo
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
        ERROR(<< "UpdateTimeStepsLength is zero.");
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

// =============
// Request Data
// ============

int LCS::RequestData(
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

    // Computations
    output->ShallowCopy(input);

    DEBUG(<< "Success");
    return 1;
}
