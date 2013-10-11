/*
 * =====================================================================================
 *
 *       Filename:  FlowMap.cxx
 *
 *    Description:  Flow Map
 *
 *        Version:  1.0
 *        Created:  01/09/2012 02:25:26 PM
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

#include <FlowMap.h>
#include <Seed.h>
#include <Interpolator.h>

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// Keys
#include <vtkInformationDoubleVectorKey.h>
#include <FilterInformation.h>


// General
#include <vtkSmartPointer.h>
#include <GeneralMath.h>

// Data
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkMultiBlockDataSet.h>

// Arrays
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

// Debug
#include <vtkInformationRequestKey.h>

// Test
#include <vtkExecutive.h>

// ======
// Macros
// ======

vtkStandardNewMacro(FlowMap);
vtkCxxRevisionMacro(FlowMap,"$Revision 1.0$");

// ===========
// Constructor
// ===========

FlowMap::FlowMap()
{
    // Member Data
    this->Dimension = 3;
    this->IntegratorMode = INTEGRATOR_MODE_USE_AdamsBashforth;
    this->IntegratorOrder = 2;
    this->IntegrationTimeStep = 1e-2;
    this->IntegrationDuration = 1;
    this->SeedReleaseGlobalTime = 0;

    // Internal Member Data
    this->IntegrationTimeStepIndex = 0;
    this->IntegrationTimeStepIndexMax = 0;
    this->IntegratorCoefficients = NULL;
    this->Tracers = NULL;

    // Pipeline
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

// ==========
// Destructor
// ==========

FlowMap::~FlowMap()
{
    // Integrator Coefficients
    if(this->IntegratorCoefficients != NULL)
    {
        delete [] this->IntegratorCoefficients;
        this->IntegratorCoefficients = NULL;
    }

    // Tracers
    if(this->Tracers != NULL)
    {
        this->Tracers->Delete();
        this->Tracers = NULL;
    }
}

// ==========
// Print Self
// ==========

void FlowMap::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ===================
// Accessors, Mutators
// ===================

void FlowMap::SetIntegratorMode(IntegratorModeType InputIntegratorMode)
{
    this->IntegratorMode = InputIntegratorMode;

    // Set default Integrator Order
    switch (this->IntegratorMode)
    {
        // Euler
        case INTEGRATOR_MODE_USE_Euler:
        {
            this->IntegratorOrder = 1;
        }

        // Runge Kutta
        case INTEGRATOR_MODE_USE_RungeKutta:
        {
            this->IntegratorOrder = 4;
        }

        // Adams bashforth
        case INTEGRATOR_MODE_USE_AdamsBashforth:
        {
            this->IntegratorOrder = 2;
        }

        // Not supported mode
        default:
        {
            ERROR(<< "Integrator mode is not supported.");
            vtkErrorMacro("Integrator mode is not supported.");
        }
    }
}

// Euler
void FlowMap::SetIntegratorModeToEuler()
{
    this->IntegratorMode = INTEGRATOR_MODE_USE_Euler;
    this->IntegratorOrder = 1;
}

// Runge Kutta
void FlowMap::SetIntegratorModeToRungeKutta()
{
    this->IntegratorMode = INTEGRATOR_MODE_USE_Euler;
    this->IntegratorOrder = 4;
}

// Adams-Bashforth
void FlowMap::SetIntegratorModeToAdamsBashforth()
{
    this->IntegratorMode = INTEGRATOR_MODE_USE_AdamsBashforth;
    this->IntegratorOrder = 2;
}

// ==========
// Get Output
// ==========

vtkStructuredGrid * FlowMap::GetOutput()
{
    return this->GetOutput(0);
}

vtkStructuredGrid * FlowMap::GetOutput(int port)
{
    return vtkStructuredGrid::SafeDownCast(this->GetOutputDataObject(port));
}

// ==========
// Set Output
// ==========

void FlowMap::SetOutput(vtkDataObject *OutputDataObject)
{
    this->SetOutput(0,OutputDataObject);
}

void FlowMap::SetOutput(int port,vtkDataObject *OutputDataObject)
{
    this->GetExecutive()->SetOutputData(port,OutputDataObject);
}

// =========
// Get Input
// =========

vtkPolyData * FlowMap::GetInput()
{
    return this->GetInput(0);
}

vtkPolyData * FlowMap::GetInput(int port)
{
    return vtkPolyData::SafeDownCast(this->GetExecutive()->GetInputData(port,0));
}

// =========
// Set Input
// =========

void FlowMap::SetInput(vtkDataObject *InputDataObject)
{
    this->SetInput(0,InputDataObject);
}

void FlowMap::SetInput(int port, vtkDataObject *InputDataObject)
{
    if(InputDataObject)
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

void FlowMap::AddInput(vtkDataObject *InputDataObject)
{
    this->AddInput(0,InputDataObject);
}

void FlowMap::AddInput(int port, vtkDataObject *InputDataObject)
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

int FlowMap::ProcessRequest(
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

    // Reuqest Information
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

int FlowMap::FillInputPortInformation(int port,vtkInformation *info)
{
    switch (port)
    {
        // Port 0 is connected to Interpolator Filter
        case 0:
        {
            info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkMultiBlockDataSet");
            DEBUG(<< "Success");
            return 1;    
        }

        // Port 1 is connected to Seed Filter
        case 1:
        {
            info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkStructuredGrid");
            DEBUG(<< "Success");
            return 1;
        }

        default:
        {
            ERROR(<< "Port number if not supported.");
            vtkErrorMacro("Port number is not supported.");
            return 0;
        }
    }
}

// ============================
// Fill Output Port Information
// ============================

int FlowMap::FillOutputPortInformation(int port, vtkInformation *info)
{
    // Output port is connected to Deformation Filter
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkStructuredGrid");
        DEBUG(<< "Success");
        return 1;
    }
    else
    {
        DEBUG(<< "Failure");
        return 0;
    }
}

// ===================
// Request Data Object
// ===================

int FlowMap::RequestDataObject(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // Ports
    for(unsigned int i=0; i<static_cast<unsigned int>(this->GetNumberOfOutputPorts()); i++)
    {
        // outputInfo
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

int FlowMap::RequestInformation(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input Info
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Check inputInfo
    if(inputInfo == NULL)
    {
        ERROR(<< "inputInfo is NULL.");
        vtkErrorMacro("inputInfo is NULL.");
        return 0;
    }

    // Output Info
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // Check outputInfo
    if(outputInfo == NULL)
    {
        ERROR(<< "outputInfo is NULL.");
        vtkErrorMacro("outputInfo is NULL.");
        return 0;
    }

    // 1- Time Steps //
    
    // Check if key exists in inputInfo
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
        vtkErrorMacro("TimeSteps is NULL");
        return 0;
    }

    // Set to Output Info
    outputInfo->Set(FilterInformation::TIME_STEPS(),TimeSteps,TimeStepsLength);

    // 2- Time Range //

    // Check if key exists in inputInfo
    if(!inputInfo->Has(FilterInformation::TIME_RANGE()))
    {
        ERROR(<< "inputInfo does not have TIME_RANGE key.");
        vtkErrorMacro("inputInfo does not have TIME");
        return 0;
    }

    double *TimeRange = inputInfo->Get(FilterInformation::TIME_RANGE());

    // Check Time Range
    if(TimeRange == NULL)
    {
        ERROR(<< "TimeRange is NULL.");
        vtkErrorMacro("TimeRange is NULL.");
        return 0;
    }

    // Set to Output Info
    outputInfo->Set(FilterInformation::TIME_RANGE(),TimeRange,2);
    
    DEBUG(<< "Success");
    return 1;
}

// =====================
// Request Update Extent
// =====================

int FlowMap::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input 0 - From Interpolator
    vtkInformation *inputInfo0 = inputVector[0]->GetInformationObject(0);
    vtkPolyData *input0 = vtkPolyData::SafeDownCast(inputInfo0->Get(vtkDataObject::DATA_OBJECT()));

    // Input 1 - From Seed
    vtkInformation *inputInfo1 = inputVector[1]->GetInformationObject(0);
    vtkStructuredGrid *input1 = vtkStructuredGrid::SafeDownCast(inputInfo1->Get(vtkDataObject::DATA_OBJECT()));

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkDataObject *output = outputInfo->Get(vtkDataObject::DATA_OBJECT());

    /*
    // We do not need TimeSteps from the output.
    // 1- TIME STEPS //

    // Check outputInfo has TIME_STEPS key
    if(!outputInfo->Has(FilterInformation::TIME_STEPS()))
    {
        ERROR(<< "outputInfo does not have TIME_STEPS key.");
        vtkErrorMacro("outputInfo does not have TIME_STEPS key.");
        return 0;
    }

    // Get TimeSteps from outputInfo
    double *TimeSteps = outputInfo->Get(FilterInformation::TIME_STEPS());
    unsigned int TimeStepsLength = outputInfo->Length(FilterInformation::TIME_STEPS());

    // Check TimeSteps
    if(TimeStepsLength < 1)
    {
        ERROR(<< "TimeSteps length is zero.");
        vtkErrorMacro("TimeSteps length is zero.");
        return 0;
    }

    if(TimeSteps == NULL)
    {
        ERROR(<< "TimeSteps is NULL.");
        vtkErrorMacro("TimeSteps is NULL.");
        return 0;
    }
    */
    

    // 2- UPDATE TIME STEPS //

    // Check outputInfo has UPDATE_TIME_STEPS key
    if(!outputInfo->Has(FilterInformation::UPDATE_TIME_STEPS()))
    {
        ERROR(<< "outoutInfo does not have UPDATE_TIME_STEPS key.");
        vtkErrorMacro("outputInfo does not have UPDATE_TIME_STEPS key.");
        return 0;
    }

    // Get Update Time Steps
    double *UpdateTimeSteps = outputInfo->Get(FilterInformation::UPDATE_TIME_STEPS());
    unsigned int UpdateTimeStepsLength = outputInfo->Length(FilterInformation::UPDATE_TIME_STEPS());

    // Check Update Time Steps
    if(UpdateTimeStepsLength < 1)
    {
        ERROR(<< "UpdateTimeSteps length is zero.");
        vtkErrorMacro("UpdateTimeSteps length is zero.");
        return 0;
    }

    if(UpdateTimeSteps == NULL)
    {
        ERROR(<< "UpdateTimeSteps is NULL.");
        vtkErrorMacro("UpdateTimeSteps is NULL.");
        return 0;
    }

    // Set UpdateTimeSteps to inputInfo
    inputInfo0->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);

    DEBUG(<< "Success");
    return 1;
}

// ========================
// Update Input Information
// ========================

// Without argument
void FlowMap::UpdateInputInformation()
{

}

// With argument
void FlowMap::UpdateInputInformation(
        vtkInformation *inputInfo,
        vtkInformation *outputInfo)
{
    // Output //

    // Check output has UPDATE_TIME_STEPS key
    if(!outputInfo->Has(FilterInformation::UPDATE_TIME_STEPS()))
    {
        ERROR(<< "outputInfo does not have UPDATE_TIME_STEPS key.");
        vtkErrorMacro("outputInfo does not have UPDATE_TIME_STEPS key.");
    }

    // Update Time Steps from output
    double *OutputUpdateTimeStep = outputInfo->Get(FilterInformation::UPDATE_TIME_STEPS());

    // Seed Release Time
    double CurrentSeedReleaseGlobalTime = OutputUpdateTimeStep[0]; // TODO for multiple release

    // Update Seed Release Time
    this->SeedReleaseGlobalTime = CurrentSeedReleaseGlobalTime;

    // Update Time Steps Length
    unsigned int UpdateTimeStepsLength = 1; // EULER for now TODO

    // Update Time Steps from input
    double UpdateTimeSteps[UpdateTimeStepsLength];

    // Set Update Time Steps
    for(unsigned int i=0; i<UpdateTimeStepsLength; i++)
    {
        UpdateTimeSteps[i] = this->GetCurrentGlobalTime(CurrentSeedReleaseGlobalTime);
    }

    // Update input Info
    inputInfo->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);
}

// ========================
// Update Input Data Object
// ========================

// Without DataObject argument
void FlowMap::UpdateInputDataObject()
{
    // Check if Traces are defined
    if(this->Tracers == NULL)
    {
        this->InitializeTracers();
    }
    else if(this->Tracers->GetNumberOfPoints() == 0)
    {
        ERROR(<< "Tracers are empty.");
        vtkErrorMacro("Tracers are empty.");
    }

    // Update input for interpolation
    this->UpdateInputDataObject(this->Tracers);
}

// With DataObject argument
void FlowMap::UpdateInputDataObject(vtkDataObject *UpdateDataObject)
{
    // Cast DataObject to PolyData
    vtkPolyData *UpdatePolyData = vtkPolyData::SafeDownCast(UpdateDataObject);

    // Get Input data object
    int Port = 0;
    int Connection = 0;
    vtkDataObject *InputDataObject = this->GetInputDataObject(Port,Connection);      ////////

    // Check input data object
    if(InputDataObject == NULL)
    {
        ERROR(<< "Input DataObject is NULL.");
        vtkErrorMacro("Input DataObject is NULL.");
    }

    // Cast Input DataObject to MultiBlock data
    vtkMultiBlockDataSet *InputMultiBlockData = vtkMultiBlockDataSet::SafeDownCast(InputDataObject);

    /// TEST ///
    
    std::cout << "IN INPUT UPDATES" << std::endl;
    if(InputMultiBlockData == NULL)
    {
        ERROR(<< "InputMultiBlockData is NULL.");
    }
    std::cout << "Number Of Blocks: " << InputMultiBlockData->GetNumberOfBlocks() << std::endl;

    ////////////

    // Set Updated DataObject to input block DataSet
    InputMultiBlockData->SetBlock(0,UpdatePolyData);

    ////////////

    std::cout << "Number Of Blocks: " << InputMultiBlockData->GetNumberOfBlocks() << std::endl;
    vtkPolyData *Test = vtkPolyData::SafeDownCast(InputMultiBlockData->GetBlock(0));
    std::cout << "Number of points: " << Test->GetNumberOfPoints() << std::endl;

    /// TEST  2 ///

    vtkInformation *InputInformation2 = this->GetInputPortInformation(0);
    if(InputInformation2 == NULL)
    {
        ERROR(<< "InputInformation2 is NULL");
    }
    vtkDataObject *InputDataObject2 = InputInformation2->Get(vtkDataObject::DATA_OBJECT());
    // vtkDataObject *InputDataObject2 = this->GetInputDataObject(Port,Connection);
    if(InputDataObject2 == NULL)
    {
        ERROR(<< "InputDataObject2 is NULL");
    }

    if(InputDataObject == InputDataObject2)
    {
        std::cout << "TWO INPUTS ARE THE SAME." << std::endl;
    }
    else
    {
        std::cout << "TWO INPUTS ARE not THE SAME." << std::endl;
    }
    vtkMultiBlockDataSet *InputMultiBlockData2 = vtkMultiBlockDataSet::SafeDownCast(InputDataObject2);
    if(InputMultiBlockData2 == NULL)
    {
        ERROR(<< "InputMultiBlockData2 is NULL.");
    }
    std::cout << "InputDataObject2: " << InputDataObject2 << std::endl;
    // InputDataObject2->Print(std::cout);
    // std::cout << "NEW: Number of Blocks2: " << InputMultiBlockData2->GetNumberOfBlocks() << std::endl;
    // vtkPolyData *InputPolyData2 = vtkPolyData::SafeDownCast(InputMultiBlockData2->GetBlock(0));
    // if(InputPolyData2 == NULL)
    // {
    //     ERROR(<< "InputPolyData2 is NULL");
    // }
    // std::cout << "New: Number of Points: " << InputPolyData2->GetNumberOfPoints() << std::endl;


    ////////////
    

    /// Test 3 ///

    vtkExecutive *Executive = this->GetExecutive();
    vtkInformationVector *InputInformationVector3 = Executive->GetInputInformation(0);
    vtkInformation *InputInformation3 = InputInformationVector3->GetInformationObject(0);
    vtkDataObject *InputDataObject3 = InputInformation3->Get(vtkDataObject::DATA_OBJECT());
    if(InputDataObject3 == InputDataObject)
    {
        std::cout << "3 and 1 are the same." << std::endl;
    }
    else
    {
        std::cout << "3 and 1 are not the same." << std::endl;
    }

}

// ===================
// Update Input Filter
// ===================

void FlowMap::UpdateInputFilter()
{
    // Update FlowMap's Input Data Object
    this->UpdateInputDataObject();
    HERE

    /// TEST ///

    vtkInformation *InputInformation = this->GetInputPortInformation(0);
    vtkDataObject *InputDataObject = InputInformation->Get(vtkDataObject::DATA_OBJECT());
    if(InputDataObject == NULL)
    {
        ERROR(<< "InputDataObject is NULL");
    }
    vtkMultiBlockDataSet *InputMultiBlockData = vtkMultiBlockDataSet::SafeDownCast(InputDataObject);
    if(InputMultiBlockData == NULL)
    {
        ERROR(<< "InputMultiBlockData is NULL.");
    }
    std::cout << "NEW: Number of Blocks: " << InputMultiBlockData->GetNumberOfBlocks() << std::endl;
    vtkPolyData *InputPolyData = vtkPolyData::SafeDownCast(InputMultiBlockData->GetBlock(0));
    if(InputPolyData == NULL)
    {
        ERROR(<< "InputPolyData is NULL");
    }
    std::cout << "New: Number of Points: " << InputPolyData->GetNumberOfPoints() << std::endl;


    ////////////

    // Update FlowMap's Input Information
    this->UpdateInputInformation();

    // Get Interpolator
    int InterpolatorPort = 0;
    int InterpolatorConnection = 0;
    vtkAlgorithm *InterpolatorAlgorithm = this->GetInputAlgorithm(InterpolatorPort,InterpolatorConnection);
    vtkDataSetAlgorithm *Interpolator = vtkDataSetAlgorithm::SafeDownCast(InterpolatorAlgorithm);

    // Modify MTime
    Interpolator->Modified();

    // Update Interpolator
    Interpolator->Update();
}

// ======================
// Set Output Data Object
// ======================

void FlowMap::SetOutputDataObject(
        vtkDataObject *SeedGridDataObject,
        vtkDataObject *TracersDataObject)
{
    // Cast Seed Grid to DataSet
    vtkDataSet *SeedGrid = vtkDataSet::SafeDownCast(TracersDataObject);

    // Cast Tracers to PolyData
    vtkPolyData *TracersPolyData = vtkPolyData::SafeDownCast(TracersDataObject);

    // Convert Tracers PolyData to Double Array
    vtkDoubleArray *TracersDataArray = vtkDoubleArray::New();
    this->ConvertPolyDataToDataArray(TracersPolyData,TracersDataArray);

    // Add DataArray to SeedGrid as vector attribute
    SeedGrid->GetPointData()->SetVectors(TracersDataArray);
    SeedGrid->GetPointData()->SetActiveVectors(TracersDataArray->GetName());

    // Get output Information
    int OutputPort = 0;
    vtkInformation *OutputInfo = this->GetOutputPortInformation(OutputPort);

    // Set output data object
    OutputInfo->Set(vtkDataObject::DATA_OBJECT(),SeedGrid);
}

// ==============================
// Convert PolyData to Data Array
// ==============================

void FlowMap::ConvertPolyDataToDataArray(
        vtkPolyData *PolyData,
        vtkDataArray *DataArray)
{
    // Default dimension
    unsigned int PointsDimension = this->Dimension;

    // Get Number of Points
    unsigned int NumberOfPoints = PolyData->GetNumberOfPoints();

    // Check PolyData
    if(NumberOfPoints < 1)
    {
        ERROR(<< "PolyData has no point.");
        vtkErrorMacro("PolyData has no point.");
    }

    // Cast DataArray to Double Array
    vtkDoubleArray *DoubleArray = vtkDoubleArray::SafeDownCast(DataArray);

    // Set DataArray
    DoubleArray->SetNumberOfComponents(PointsDimension);
    DoubleArray->SetNumberOfTuples(NumberOfPoints);
    DoubleArray->SetName("Tracers");

    // Get Points form PolyData
    vtkPoints *Points = PolyData->GetPoints();

    // Set Progress
    this->ProgressReset();
    this->SetProgressMessage("Save output");

    // Convert points to doubles
    for(unsigned int PointIterator = 0; PointIterator < NumberOfPoints; PointIterator++)
    {
        // Get a point as double
        double *DoublePoint = Points->GetPoint(PointIterator);

        // Dimension consideratrion
        if(PointsDimension == 2)
        {
            // Set Double Array
            DoubleArray->SetTuple2(PointIterator,DoublePoint[0],DoublePoint[1]);
        }
        else if(PointsDimension == 3)
        {
            // Set Double Array
            DoubleArray->SetTupleValue(PointIterator,DoublePoint);
        }
        else
        {
            ERROR(<< "Dimension is not supported.");
            vtkErrorMacro("Dimension is not supported.");
            break;
        }

        // Update progress
        this->ProgressUpdate(PointIterator,NumberOfPoints);
    }

    // Reset progress
    this->ProgressReset();
}

// ================
// Information Test
// ================

void FlowMap::InformationTest(vtkInformationVector **inputVector)
{
    // Form inputVector
    vtkInformation *InputInformation1 = inputVector[0]->GetInformationObject(0);
    vtkDataObject *InputDataObject1 = InputInformation1->Get(vtkDataObject::DATA_OBJECT());
    vtkMultiBlockDataSet *InputMultiBlock1 = vtkMultiBlockDataSet::SafeDownCast(InputDataObject1);

    if(!InputDataObject1)
    {
        ERROR(<< "InputDataObject1 is null.");
    }
    if(!InputMultiBlock1)
    {
        ERROR(<< "InputMultiBlock1 is NULL.");
    }

    // From Executive
    vtkExecutive *Executive = this->GetExecutive();
    vtkInformationVector *InputVector2 = Executive->GetInputInformation(0);
    vtkInformation *InputInformation2 = InputVector2->GetInformationObject(0);
    vtkDataObject *InputDataObject2 = InputInformation2->Get(vtkDataObject::DATA_OBJECT());
    vtkMultiBlockDataSet *InputMultiBlock2 = vtkMultiBlockDataSet::SafeDownCast(InputDataObject2);

    if(!InputDataObject2)
    {
        ERROR(<< "InputDataObject2 is NULL.")
    }
    if(!InputMultiBlock2)
    {
        ERROR(<< "InputMultiBlock2 is NULL.")
    }

    // From Get Input Data
    vtkDataObject *InputDataObject3 = this->GetInputDataObject(0,0);
    vtkMultiBlockDataSet *InputMultiBlock3 = vtkMultiBlockDataSet::SafeDownCast(InputDataObject3);

    if(!InputDataObject3)
    {
        ERROR(<< "InputDataObject3 is NULL.");
    }
    if(!InputMultiBlock3)
    {
        ERROR(<< "InputMultiBlock3 is NULL.");
    }

    // Comparision
    std::cout << "DO1: " << InputDataObject1 << ", DO2: " << InputDataObject2 << ", DO3: " << InputDataObject3 << std::endl;
    std::cout << "MB1: " << InputMultiBlock1 << ", MB2: " << InputMultiBlock2 << ", MB3: " << InputMultiBlock3 << std::endl;
}

// ============
// Request Data
// ============

int FlowMap::RequestData(
        vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input 0 - From Interpolator
    vtkInformation *inputInfo0 = inputVector[0]->GetInformationObject(0);
    vtkPolyData *TracersData = vtkPolyData::SafeDownCast(inputInfo0->Get(vtkDataObject::DATA_OBJECT()));

    // Input 1 - From Seed
    vtkInformation *inputInfo1 = inputVector[1]->GetInformationObject(0);
    vtkStructuredGrid *SeedGrid = vtkStructuredGrid::SafeDownCast(inputInfo1->Get(vtkDataObject::DATA_OBJECT()));

    if(SeedGrid == NULL)
    {
        ERROR(<< "Seed Grid is NULL");
        vtkErrorMacro("Seed Grid is NULL");
        return 0;
    }
    else if(SeedGrid->GetNumberOfPoints() == 0)
    {
        ERROR(<< "Number of Seed Grid points is zero.");
        vtkErrorMacro("Number of Seed Grid points is zero.");
        return 0;
    }

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkStructuredGrid *output = vtkStructuredGrid::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    /*
    if(this->IntegrationTimeStepIndex == 0)
    {
        // 1- Initialization
        this->InitializeTimeIndices(inputInfo0);
        this->InitializeIntegratorCoefficients();
        this->InitializeTracers(input1);
    }
    else
    {
        // 2- Processing
        this->AdvectTracers(input0);

        // Count number of processes
        this->IntegrationTimeStepIndex++;
    }
    
    // check for Process termination
    if(this->IntegrationTimeStepIndex < this->IntegrationTimeStepIndexMax)
    {
        // Continue Processing
        // request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(),1);

        // // TEST Interpolator //
        // vtkAlgorithm *InterpolatorAlgorithm = this->GetInputAlgorithm(0,0);
        // Interpolator *InterpolatorFilter = Interpolator::SafeDownCast(InterpolatorAlgorithm);
        // std::cout << "Name: " << InterpolatorFilter->GetClassName() << std::endl;
        // vtkExecutive *InterpolatorExecutive = InterpolatorFilter->GetExecutive();
        // vtkInformationVector *InterpolatorOutputInfoVector = InterpolatorExecutive->GetOutputInformation();
        // vtkInformation *InterpolatorOutputInfo = InterpolatorOutputInfoVector->GetInformationObject(0);
        // vtkDataObject *InterpolatorOutputDataObject = InterpolatorOutputInfo->Get(vtkDataObject::DATA_OBJECT());
        // vtkMultiBlockDataSet *MultiBlockDataSet = vtkMultiBlockDataSet::SafeDownCast(InterpolatorOutputDataObject);
        // vtkPolyData *InterpolatorTracer = vtkPolyData::SafeDownCast(MultiBlockDataSet->GetBlock(0));
        // std::cout << "Numer of points: " << InterpolatorTracer->GetNumberOfPoints() << std::endl;

        // // Streaming
        // vtkStreamingDemandDrivenPipeline *InterpolatorStreamingExecutive = vtkStreamingDemandDrivenPipeline::SafeDownCast(InterpolatorExecutive);

        // InterpolatorFilter->UpdateInformation();
        // InterpolatorFilter->UpdateWholeExtent();
        // std::cout << "Interpolator before modified" << std::endl;
        // std::cout << "\tInterpolator MTime: " << InterpolatorFilter->GetMTime() << std::endl;
        // std::cout << "\tInterpolator Executive MTime: " << InterpolatorExecutive->GetMTime() << std::endl;
        // std::cout << "\tInterpolatoe Pipeline MTime: " << InterpolatorStreamingExecutive->GetPipelineMTime() << std::endl;

        // // Should not change MTime
        // InterpolatorFilter->Modified();

        // std::cout << "Interpolator after modified." << std::endl;
        // std::cout << "\tInterpolator MTime: " << InterpolatorFilter->GetMTime() << std::endl;
        // std::cout << "\tInterpolator Executive MTime: " << InterpolatorExecutive->GetMTime() << std::endl;
        // std::cout << "\tInterpolatoe Pipeline MTime: " << InterpolatorStreamingExecutive->GetPipelineMTime() << std::endl;

        // // InterpolatorStreamingExecutive->Update();
        // InterpolatorFilter->Update(0);

        // std::cout << "Interpolator after update" << std::endl;
        // std::cout << "\tInterpolator MTime: " << InterpolatorFilter->GetMTime() << std::endl;
        // std::cout << "\tInterpolator Executive MTime: " << InterpolatorExecutive->GetMTime() << std::endl;
        // std::cout << "\tInterpolatoe Pipeline MTime: " << InterpolatorStreamingExecutive->GetPipelineMTime() << std::endl;

        // // Test Seed //
        // vtkAlgorithm *SeedAlgorithm = this->GetInputAlgorithm(1,0);
        // Seed *SeedFilter = Seed::SafeDownCast(SeedAlgorithm);
        // // Seed *SeedFilter = SeedAlgorithm;
        // std::cout << "Name: " << SeedFilter->GetClassName() << std::endl;
        // vtkExecutive *SeedExecutive = SeedFilter->GetExecutive();
        // vtkInformation *SeedInputInfo = SeedExecutive->GetInputInformation(0,0);
        // vtkInformation *SeedOutputInfo = SeedExecutive->GetOutputInformation(0);

        // // Streaming
        // vtkStreamingDemandDrivenPipeline *SeedStreamingExecutive = vtkStreamingDemandDrivenPipeline::SafeDownCast(SeedExecutive);

        // // Make change in input information
        // double *SeedUpdateTimeSteps = SeedInputInfo->Get(FilterInformation::UPDATE_TIME_STEPS());
        // std::cout << "Seed update time steps: " << SeedUpdateTimeSteps[0] << std::endl;
        // SeedUpdateTimeSteps[0] = 2;

        // std::cout << "Seed before update: " << std::endl;
        // std::cout << "Seed MTime: " << SeedFilter->GetMTime() << std::endl;
        // std::cout << "Seed Executive MTime: " << SeedFilter->GetExecutive()->GetMTime() << std::endl;
        // std::cout << "\tSeed Pipeline MTime: " << SeedStreamingExecutive->GetPipelineMTime() << std::endl;

        // std::cout << "Seed before update" << std::endl;
        // // SeedFilter->UpdateWholeExtent();
        // // SeedFilter->UpdateInformation();

        // // Better
        // SeedStreamingExecutive->Update();
        // // SeedFilter->Update();
        // //
        // std::cout << "Seed after update" << std::endl;
        // std::cout << "Seed MTime: " << SeedFilter->GetMTime() << std::endl;
        // std::cout << "Seed Executive MTime: " << SeedFilter->GetExecutive()->GetMTime() << std::endl;
        // std::cout << "\tSeed Pipeline MTime: " << SeedStreamingExecutive->GetPipelineMTime() << std::endl;
        this->UpdateInputFilter();



    }
    else
    {
        // 3- Post Processing
        this->PostProcessingTracers(input0,input1,output);

        // Terminate Run
        // request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());

        // Reset Integration Counter
        this->IntegrationTimeStepIndex = 0;
    }
    */

    this->InformationTest(inputVector);

    /*
    // Initialization //

    // Initialize time indices
    this->InitializeTimeIndices(inputInfo0);

    // Initialize Integrator Coefficients
    this->InitializeIntegratorCoefficients();

    // Initialize Integrator Coefficients
    this->InitializeIntegratorCoefficients();

    // Initialize Tracers
    this->InitializeTracers();

    // Advecting Tracers //

    // Set Progress
    this->ProgressReset();
    this->SetProgressMessage("Integration steps");

    // Time step index iteration
    while(this->IntegrationTimeStepIndex < IntegrationTimeStepIndexMax)
    {
        // Acquire Tracers data
        HERE

        /// TEST ///

        vtkInformation *InputInformation = this->GetInputPortInformation(0);
        vtkDataObject *InputDataObject = InputInformation->Get(vtkDataObject::DATA_OBJECT());
        HERE
        vtkMultiBlockDataSet *InputBlocks = vtkMultiBlockDataSet::SafeDownCast(InputDataObject);
        HERE
        if(InputBlocks == NULL)
        {
            ERROR(<< "InputBlocks are NULL");
        }

        // std::cout << "Number of Blocks: " << InputBlocks->GetNumberOfBlocks() << std::endl;
        HERE

        ////////////

        this->UpdateInputFilter();
        HERE




        // Integrate
        this->AdvectTracers(TracersData);

        // Update Progress
        this->ProgressUpdate(this->IntegrationTimeStepIndex,IntegrationTimeStepIndexMax);
    }

    // Reset progress
    this->ProgressReset();

    // Set Seed Grid to output
    this->SetOutputDataObject(SeedGrid,TracersData);
    */

    DEBUG(<< "Success");
    return 1;
}

// ======================
// Integrator Mode String
// ======================

const char * FlowMap::IntegratorModeString(IntegratorModeType IntegratorMode)
{
    switch(IntegratorMode)
    {
        // Euler
        case INTEGRATOR_MODE_USE_Euler:
        {
            return "INTEGRATOR_MODE_USE_Euler";
        }

        // Runge Kutta 4
        case INTEGRATOR_MODE_USE_RungeKutta:
        {
            return "INTEGRATOR_MODE_USE_RungeKutta";
        }

        // Adams Bashforth
        case INTEGRATOR_MODE_USE_AdamsBashforth:
        {
            return "INTEGRATOR_MODE_USE_Adamsbashforth";
        }

        // Not supported mode
        default:
        {
            ERROR(<< "Integrator mode is not supported.");
            vtkErrorMacro("Integrator mode is not supported.");
            return NULL;
        }
    }
}

// ====================
// Required Time Frames
// ====================

unsigned int FlowMap::RequiredTimeFrames(IntegratorModeType ModeType)
{
    switch (ModeType)
    {
        // Euler
        case INTEGRATOR_MODE_USE_Euler:
        {
            return 1;
        }

        // Runge-Kutta 4
        case INTEGRATOR_MODE_USE_RungeKutta:
        {
            // To be eddited
            return this->IntegratorOrder;
        }

        // Adams-Bashforth
        case INTEGRATOR_MODE_USE_AdamsBashforth:
        {
            // To be editted
            return this->IntegratorOrder;
        }

        // Not supported Mode
        default:
        {
            ERROR(<< "Integrator mode not supported.");
            vtkErrorMacro("Integrator mode not supported.");
            return 1;
        }
    }
}

// =======================
// Initialize Time Indices
// =======================

void FlowMap::InitializeTimeIndices(vtkInformation *inputInfo0)
{
    // Get Data Time Range
    double *DataTimeRange = inputInfo0->Get(FilterInformation::TIME_RANGE());
    double DataTimeRangeInterval = DataTimeRange[1] - DataTimeRange[0];

    // Check Integration Duration
    if(this->IntegrationDuration > DataTimeRangeInterval)
    {
        ERROR(<< "Integration Duration is bigger than Data Time Rang Interval.");
        vtkErrorMacro("Intergation Duration is bigger than Data Time Range Interval.");
    }

    // Set Integration Time Step Index Max
    double DurationToTimeStepRatio = (double) this->IntegrationDuration / (double) this->IntegrationTimeStep;

    // Check ratio is integer
    if(floor(DurationToTimeStepRatio) != DurationToTimeStepRatio)
    {
        WARNING(<< "IntegrationDuration is not multiple of Integration Time Step.");
        DurationToTimeStepRatio = round(DurationToTimeStepRatio);
    }

    this->IntegrationTimeStepIndexMax = static_cast<unsigned long int>(DurationToTimeStepRatio);
}

// ==================================
// Initialize Integrator Coefficients
// ==================================

void FlowMap::InitializeIntegratorCoefficients()
{
    // Initialize array
    this->IntegratorCoefficients = new double[this->IntegratorOrder];

    // Use specific integrator method
    switch(this->IntegratorMode)
    {
        // 1- Runge-Kutta //
        case INTEGRATOR_MODE_USE_RungeKutta:
        {
            // TODO
            break;
        }

        // 2- Adams-Bashforth //
        case INTEGRATOR_MODE_USE_AdamsBashforth:
        {
            // Set Order for Adams-Bashforth
            unsigned int order = this->IntegratorOrder;

            // Allocate Rational Coefficients
            RationalNumber *RationalCoefficients = new RationalNumber[this->IntegratorOrder];

            // Compute Rational Coefficients
            this->AdamsBashforthCoefficients(RationalCoefficients,order);

            // Convert Rational Numbers to double type
            for(unsigned int i=0; i<this->IntegratorOrder; i++)
            {
                this->IntegratorCoefficients[i] = (*(RationalCoefficients+i)).Evaluate();
            }

            // Delete Rational Coefficients
            delete [] RationalCoefficients;

            break;
        }

        // 3- Not supported Integrator mode //
        default:
        {
            DEBUG(<< "Integrator mode is not suppoeted.");
            vtkErrorMacro("Integrator mode is not supported.");
        }
    }
}

// ============================
// Adams-Bashforth Coefficients
// ============================

void FlowMap::AdamsBashforthCoefficients(RationalNumber *RationalCoefficients,unsigned int order)
{
    // Iterate over each coefficient
    for(unsigned int j=0; j<order; j++)
    {
        // Initialize Polynomial P
        Polynomial P(1);

        // Create Polynomial P
        for(unsigned int i=0; i<order; i++)
        {
            // Define Polynomial Q
            Polynomial Q;

            // If i==i, Q = 1
            if(i == j)
            {
                Q.SetDegree(0);
                Q.SetCoefficient(0,1);
            }
            // If i != j, Q = ("x"+i)
            else
            {
                Q.SetDegree(1);
                Q.SetCoefficient(0,i);
                Q.SetCoefficient(1,1);
            }

            // Successive Multiplication
            P *= Q;
        }

        // Integrate P
        Polynomial IntP = P.Integrate();
        RationalNumber IntPValue = IntP.Evaluate(RationalNumber(1))-IntP.Evaluate(RationalNumber(0));

        // Fraction
        RationalNumber FractionNumerator = RationalNumber(-1)^static_cast<int>(j);
        RationalNumber FractionDenominator1 = Combinatorics::Factorial(j);
        RationalNumber FractionDenominator2 = Combinatorics::Factorial(order-j-1);
        RationalNumber Fraction = FractionNumerator / (FractionDenominator1 * FractionDenominator2);

        // Coefficient
        RationalCoefficients[order-j-1] = Fraction * IntPValue;
    }
}

// ========================
// Runge-Kutta Coefficinets
// ========================

void  FlowMap::RungeKuttaCoefficients(RationalNumber *RationalCoefficients, unsigned int order)
{
    // TODO
}

// ==================
// Initialize Tracers
// ==================

// Without argument
void FlowMap::InitializeTracers()
{
    // Seed Input Port/Connection
    unsigned int SeedFilterPort = 1;
    unsigned int SeedFilterConnection = 0;

    // Filter's Input Information
    vtkInformation *InputInformation = this->GetInputInformation(SeedFilterPort,SeedFilterConnection);

    // Filter's Input DataObject
    vtkDataObject *SeedDataObject = InputInformation->Get(vtkDataObject::DATA_OBJECT());

    // Seed's StructuredGrid output
    vtkStructuredGrid *SeedGrid = vtkStructuredGrid::SafeDownCast(SeedDataObject);

    // Check Seed Grid
    if(SeedGrid == NULL || SeedGrid->GetNumberOfPoints() < 1)
    {
        // Prepare Seed Grid //

        // Get Seed Filter
        vtkAlgorithm *SeedAlgorithm = this->GetInputAlgorithm(SeedFilterPort,SeedFilterConnection);

        // Cast to Seed DaataSet Algorithm
        vtkDataSetAlgorithm *SeedFilter = vtkDataSetAlgorithm::SafeDownCast(SeedAlgorithm);

        // Modify Seed MTime
        SeedFilter->Modified();

        // Update Seed
        SeedFilter->Update();
    }

    // Check Seed Grid again
    if(SeedGrid == NULL || SeedGrid->GetNumberOfPoints() < 1)
    {
        ERROR(<< "SeedGrid is empty or has no points.");
        vtkErrorMacro("SeedGrid is empty or has no points.");
    }

    // Initialize Tracers
    this->InitializeTracers(SeedGrid);
}

// With argument
void FlowMap::InitializeTracers(vtkStructuredGrid *SeedGrid)
{
    // Initialize Tracers
    this->Tracers = vtkPolyData::New();
    
    // SeedGrid Points
    vtkPoints *SeedPoints = SeedGrid->GetPoints();

    // Check SeedPoints
    if(SeedPoints == NULL)
    {
        ERROR(<< "SeedPoints is NULL.");
        vtkErrorMacro("SeedPoints is NULL.");
    }
    else if(SeedPoints->GetNumberOfPoints() < 1)
    {
        ERROR(<< "SeedPoints has no points.");
        vtkErrorMacro("SeedPoints has no points.");
    }

    // Tracer Points
    vtkPoints *TracerPoints = vtkPoints::New();
    if(TracerPoints == NULL)
    {
        ERROR(<< "TracerPoints is NULL.");
    }

    // Copy Points from Seed to Tracers
    TracerPoints->DeepCopy(SeedPoints);

    // Set Points of Tracer
    this->Tracers->SetPoints(TracerPoints);
}

// =======================
// Get Current Global Time
// =======================

double FlowMap::GetCurrentGlobalTime(double CurrentSeedReleaseGlobalTime) const
{
    double CurrentGlobalTime = CurrentSeedReleaseGlobalTime +
        this->IntegrationTimeStepIndex * IntegrationTimeStep;

    return CurrentGlobalTime;
}

// ==============
// Advect Tracers
// ==============

void FlowMap::AdvectTracers(vtkPolyData *TracersData)
{
    // Check Tracers
    if(this->Tracers == NULL)
    {
        ERROR(<< "Tracers are NULL.");
        vtkErrorMacro("Tracers are NULL.");
    }

    // Check TracerData
    if(TracersData == NULL)
    {
        ERROR(<< "Tracers Data is NULL.");
        vtkErrorMacro("Tracers Data is NULL.");
    }

    // Number of Tracers
    unsigned int NumberOfTracers = this->Tracers->GetNumberOfPoints();

    HERE

    // Tracer Points
    vtkPoints *TracerPoints = this->Tracers->GetPoints();

    // Velocities
    vtkDoubleArray *Velocities = vtkDoubleArray::SafeDownCast(TracersData->GetPointData()->GetVectors("velocity"));
    std::cout << Velocities->GetNumberOfTuples() << std::endl;

    DISPLAY(Velocities->GetNumberOfComponents());
    DISPLAY(Velocities->GetNumberOfTuples());

    // Declare Point and Velocity array for each point
    double Point[this->Dimension];
    double Velocity[this->Dimension];

    // Set Progress
    this->ProgressReset();
    this->SetProgressMessage("Advect Tracers");

    // Iterate over points
    std::cout << "In flowmap::Advect" << std::endl;
    for(unsigned int i=0; i<NumberOfTracers; i++)
    {
        // Get a Point
        TracerPoints->GetPoint(i,Point);

        // Get Velocity
        Velocities->GetTuple(i,Velocity);

        // Forward Euler Method
        for(unsigned int j=0; j<this->Dimension; j++)
        {
            Point[j] += Velocity[j] * this->IntegrationTimeStep;
        }

        // Update Progress
        this->ProgressUpdate(i,NumberOfTracers);
    }

    // Reset Progress
    this->ProgressReset();
}

// =======================
// Post Processing Tracers
// =======================

void FlowMap::PostProcessingTracers(
        vtkPolyData *input0,
        vtkStructuredGrid *input1,
        vtkStructuredGrid *output)
{
    //
}
