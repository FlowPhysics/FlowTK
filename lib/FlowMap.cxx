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
    // Check if Traces as defined
    if(this->Tracers == NULL)
    {
        // this->InitializeTracers();
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
    vtkDataObject *InputDataObject = this->GetInputDataObject(Port,Connection);

    // Check input data object
    if(InputDataObject == NULL)
    {
        ERROR(<< "Input DataObject is NULL.");
        vtkErrorMacro("Input DataObject is NULL.");
    }

    // Cast Input DataObject to MultiBlock data
    vtkMultiBlockDataSet *InputMultiBlockData = vtkMultiBlockDataSet::SafeDownCast(InputDataObject);

    // Set Updated DataObject to input block DataSet
    InputMultiBlockData->SetBlock(0,UpdatePolyData);
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
    vtkPolyData *input0 = vtkPolyData::SafeDownCast(inputInfo0->Get(vtkDataObject::DATA_OBJECT()));

    // Input 1 - From Seed
    vtkInformation *inputInfo1 = inputVector[1]->GetInformationObject(0);
    vtkStructuredGrid *input1 = vtkStructuredGrid::SafeDownCast(inputInfo1->Get(vtkDataObject::DATA_OBJECT()));

    if(input1 == NULL)
    {
        ERROR(<< "input1 is NULL");
        vtkErrorMacro("input1 is NULL");
        return 0;
    }
    else if(input1->GetNumberOfPoints() == 0)
    {
        ERROR(<< "Number of input1 points is zero.");
        vtkErrorMacro("Number of input1 points is zero.");
        return 0;
    }

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkStructuredGrid *output = vtkStructuredGrid::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Computation //

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
        // Update input of port 0
        this->UpdateInputDataObject();
        this->UpdateInputInformation(inputInfo0,outputInfo);
        
        // Continue Processing
        request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(),1);
    }
    else
    {
        // 3- Post Processing
        this->PostProcessingTracers(input0,input1,output);

        // Terminate Run
        request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());

        // Reset Integration Counter
        this->IntegrationTimeStepIndex = 0;
    }

    // Update output
    output->ShallowCopy(input1);

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

double FlowMap::GetCurrentGlobalTime(double CurrentSeedReleaseGlobalTime)
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
