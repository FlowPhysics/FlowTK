/*
 * =====================================================================================
 *
 *       Filename:  Interpolator.cxx
 *
 *    Description:  Interpolates for requested time steps
 *
 *        Version:  1.0
 *        Created:  08/31/2012 05:04:03 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of Calfornia, Berkeley
 *
 * =====================================================================================
 */

// =======
// Headers
// =======

// For Pipeline
#include <Interpolator.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// General
#include <vtkSmartPointer.h>

// Key
#include <vtkInformationDoubleVectorKey.h>
#include <FilterInformation.h>

// Data
#include <vtkCompositeDataIterator.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>

// Elements
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkGenericCell.h>

// Arrays
#include <vtkIdTypeArray.h>
#include <vtkDoubleArray.h>

// For DEBUG
#include <vtkInformationRequestKey.h>

// ======
// Macros
// ======

#define EPSILON 1e-4 // for cell search
vtkStandardNewMacro(Interpolator);
vtkCxxRevisionMacro(Interpolator,"$Revision 1.0$");

// ===========
// Constructor
// ===========

Interpolator::Interpolator()
{
    // Member Data
    this->SnapToTimeStepTolerance = 1e-4;

    // Internal Member data
    this->Dimension = 3;
    this->GridAdjacencies = NULL;

    // Pipeline
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

// ==========
// Destructor
// ==========

Interpolator::~Interpolator()
{
    // Grid Adjacencies
    if(this->GridAdjacencies != NULL)
    {
        this->GridAdjacencies->Delete();
        this->GridAdjacencies = NULL;
    }
}

// ==========
// Print Self
// ==========

void Interpolator::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ===================
// Accessors, Mutators
// ===================


// ==========
// Get Output
// ==========

vtkUnstructuredGrid * Interpolator::GetOutput()
{
    return this->GetOutput(0);
}

vtkUnstructuredGrid * Interpolator::GetOutput(int port)
{
    return vtkUnstructuredGrid::SafeDownCast(this->GetOutputDataObject(port));
}

// ==========
// Set Output
// ==========

void Interpolator::SetOutput(vtkDataObject *OutputDataObject)
{
    this->SetOutput(0,OutputDataObject);
}

void Interpolator::SetOutput(int port, vtkDataObject *OutputDataObject)
{
    this->GetExecutive()->SetOutputData(port,OutputDataObject);
}

// =========
// Get Input
// =========

vtkDataSet * Interpolator::GetInput()
{
    return this->GetInput(0);
}

vtkDataSet * Interpolator::GetInput(int port)
{
    return vtkDataSet::SafeDownCast(this->GetExecutive()->GetInputData(port,0));
}

// =========
// Set Input
// =========

void Interpolator::SetInput(vtkDataObject *InputDataObject)
{
    this->SetInput(0,InputDataObject);
}

void Interpolator::SetInput(int port, vtkDataObject *InputDataObject)
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

void Interpolator::AddInput(vtkDataObject *InputDataObject)
{
    this->AddInput(0,InputDataObject);
}

void Interpolator::AddInput(int port, vtkDataObject *InputDataObject)
{
    if(InputDataObject != NULL)
    {
        #if VTK_MAJOR_VERSION <= 5
        this->AddInputConnection(port,InputDataObject->GetProducerPort());
        #else
        this->AddInputDataObject(port,InputDataObject);
        #endif
    }
}

// ===============
// Process Request
// ===============

int Interpolator::ProcessRequest(
        vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    DEBUG(<< request->GetRequest()->GetName());
    // request->Print(std::cout);

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

    // Reuqest Update Extent
    if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
        return this->RequestUpdateExtent(request,inputVector,outputVector);
    }

    // Request Data
    if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
        return this->RequestData(request,inputVector,outputVector);
    }

    // Request Data Not Generated
    if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_NOT_GENERATED()))
    {
        unsigned int OutputPort = request->Get(vtkDemandDrivenPipeline::FROM_OUTPUT_PORT());
        vtkInformation *OutputInfo = outputVector->GetInformationObject(OutputPort);
        OutputInfo->Set(vtkDemandDrivenPipeline::DATA_NOT_GENERATED(),1);
    }

    // Otherwise use supperclass
    return this->Superclass::ProcessRequest(request,inputVector,outputVector);
}

// ===========================
// Fill Input Port Information
// ===========================

int Interpolator::FillInputPortInformation(int port,vtkInformation *info)
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

int Interpolator::FillOutputPortInformation(int port, vtkInformation *info)
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

int Interpolator::RequestDataObject(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // Prots
    for(unsigned int i=0; i<static_cast<unsigned int>(this->GetNumberOfOutputPorts()); i++)
    {
        HERE
        this->Test();

        // Output
        vtkInformation *outputInfo = outputVector->GetInformationObject(i);
        vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

        // Create new instance
        if(output == NULL)
        {
            output = vtkMultiBlockDataSet::New();
            outputInfo->Set(vtkDataObject::DATA_OBJECT(),output);
            output->FastDelete();
            this->GetOutputPortInformation(i)->Set(vtkDataObject::DATA_EXTENT_TYPE(),output->GetExtentType());
        }

        HERE
        this->Test();
    }

    DEBUG(<< "Success");
    return 1;
}

// ===================
// Request Information
// ===================

int Interpolator::RequestInformation(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    HERE
    this->Test();

    // Input 
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // 1- TIME STEPS //
    
    // Check if Key exists in inputInfo
    if(!inputInfo->Has(FilterInformation::TIME_STEPS()))
    {
        ERROR(<< "inputInfo does not have TIME_STEPS key.");
        vtkErrorMacro("inputInfo does not have TIME_STEPS key.");
        return 1;
    }

    // Get Time Steps from Input
    unsigned int TimeStepsLength = inputInfo->Length(FilterInformation::TIME_STEPS());
    double *TimeSteps = inputInfo->Get(FilterInformation::TIME_STEPS());

    // Check Time Steps
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

    // Set Time Steps to Output
    outputInfo->Set(FilterInformation::TIME_STEPS(),TimeSteps,TimeStepsLength);

    // 2- TIME RANGE //

    // Check if TimeRange key exists in inputInfo
    if(!inputInfo->Has(FilterInformation::TIME_RANGE()))
    {
        ERROR(<< "inputInfo does not have TimeRange key.");
        vtkErrorMacro("inputInfo does not have TimeRange key.");
        return 0;
    }

    // Get Time Range from Input
    double *TimeRange = inputInfo->Get(FilterInformation::TIME_RANGE());

    // Check TimeRange
    if(TimeRange == NULL)
    {
        ERROR(<< "TimeRange is NULL.");
        vtkErrorMacro("TimeRange is NULL.");
        return 0;
    }

    // Set Time Range to Output
    outputInfo->Set(FilterInformation::TIME_RANGE(),TimeRange,2);

    HERE
    this->Test();

    DEBUG(<< "Success");
    return 1;
}

// =====================
// Request Update Extent
// =====================

int Interpolator::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input Information
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Check inputInfo
    if(inputInfo == NULL)
    {
        ERROR(<< "inputInfo is NULL.");
        vtkErrorMacro("inputInfo is NULL.");
        return 0;
    }

    // Output Information
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // Check outputInfo
    if(outputInfo == NULL)
    {
        ERROR(<< "outputInfo is NULL.");
        vtkErrorMacro("outputInfo is NULL.");
        return 0;
    }

    // 1- DATA TIME STEPS from Input //
    
    // Check inputInfo has DATA_TIME_STEPS key
    if(!inputInfo->Has(FilterInformation::DATA_TIME_STEPS()))
    {
        ERROR(<< "inputInfo does not have DATA_TIME_STEPS key.");
        vtkErrorMacro("inputInfo does not have DATA_TIME_STEPS key.");
        return 0;
    }

    // Get Time Steps from Input
    unsigned int InputDataTimeStepsLength = inputInfo->Length(FilterInformation::DATA_TIME_STEPS());
    double *InputDataTimeSteps = inputInfo->Get(FilterInformation::DATA_TIME_STEPS());

    // Check Input Data TimeSteps
    if(InputDataTimeStepsLength < 1)
    {
        ERROR(<< "InputDataTimeSteps length is zero.");
        vtkErrorMacro("InputDataTimeSteps length is zero.");
        return 0;
    }

    if(InputDataTimeSteps == NULL)
    {
        ERROR(<< "InputDataTimeSteps is NULL.");
        vtkErrorMacro("InputDataTimeSteps is NULL.");
        return 0;
    }

    DISPLAY(InputDataTimeSteps,InputDataTimeStepsLength);

    // 2- UPDATE TIME STEPS from Output //
    
    // Check outputInfo has UPDATE_TIME_STEPS key
    if(!outputInfo->Has(FilterInformation::UPDATE_TIME_STEPS()))
    {
        ERROR(<< "outputInfo does not have UPDATE_TIME_STEPS key.");
        vtkErrorMacro("outputInfo does not have UPDATE_TIME_STEPS key.");
        return 0;
    }

    // Get Update Time Step from Output
    unsigned int OutputUpdateTimeStepsLength = outputInfo->Length(FilterInformation::UPDATE_TIME_STEPS());
    double *OutputUpdateTimeSteps = outputInfo->Get(FilterInformation::UPDATE_TIME_STEPS());

    // Check UpdateTimeSteps
    if(OutputUpdateTimeStepsLength < 1)
    {
        WARNING(<< "No update requested.");
        vtkErrorMacro("No update requested.");
        return 0;
    }

    DISPLAY(OutputUpdateTimeSteps,OutputUpdateTimeStepsLength);

    // 3- UPDATE TIME STEPS FOR INPUT //

    // Requested data TimeStep Intervals
    unsigned int RequestedDataTimeIntervals = OutputUpdateTimeSteps;
    unsigned int *RequestedDataTimeStepIndices[RequestedDataTimeIntervalsLength];
    unsigned int *InterpolationRequested[RequestedDataTimeStepsIntervalLength];

    this->FindRequestedDataTimeIntervals(
            InputDataTimeSteps,
            InputDataTimeStepsLength,
            OutputUpdateTimeSteps,
            OutputUpdateTimeStepsLength,
            RequestedDataTimeIntervals,          // Output
            InterpolationRequested);             // Output

    this->FindRequestedDataTimeStepsIndices(
            RequestedDataTimeIntervals,
            RequestedDataTimeIntervalsLength,
            RequestedDataTimeStepIndices,         //Output
            RequestedDataTimeStepIndicesLength);  // Output
    

//    DISPLAY(SqueezedRequestIndex,NumberOfSqueezedIndices);

    // Check sorted
    if(NumberOfSqueezedIndices>1)
    {
        for(unsigned int i=0; i<NumberOfSqueezedIndices-1; i++)
        {
            if(SqueezedRequestIndex[i]>SqueezedRequestIndex[i+1])
            {
                vtkErrorMacro("Interpolator >> RequestUpdateExtent >> SqueezedRequestIndex are not sorted.");
                return 0;
            }
        }
    }

    // Convert Indices to Time Steps
    double InputUpdateTimeSteps[NumberOfSqueezedIndices];

    for(unsigned int i=0; i<NumberOfSqueezedIndices; i++)
    {
        InputUpdateTimeSteps[i] = InputTimeSteps[SqueezedRequestIndex[i]];
    }

    DISPLAY(InputUpdateTimeSteps,NumberOfSqueezedIndices)

    // Set Input Update Time Steps to Input
    inputInfo->Set(FilterInformation::UPDATE_TIME_STEPS(),InputUpdateTimeSteps,NumberOfSqueezedIndices);

    return 1;
}

// =======================================
// Find Requested Data Time Step Intervals
// =======================================

bool Interpolator::FindRequestedDataTimeIntervals(
        double *InputDataTimeSteps,
        unsigned int InputDataTimeStepsLength,
        double *OutputUpdateTimeSteps,
        unsigned int OutputUpdateTimeStepsLength,
        unsigned int *RequestedDataTimeIntervals,       // Output
        bool *InterpolationRequested)                   // Output
{
    // Status of finding intervals
    bool IntervalFound;

    // Loop over OutputUpdateTimeSteps
    for(unsigned int OutputUpdateTimeStepsIterator = 0;
        OutputUpdateTimeStepsIterator < OutputUpdateTimeStepsLength;
        OutputUpdateTimeStepsIterator++, IntervalFound = false)
    {
        // Check if UpdateTimeSteps are in the Time Range
        if(OutputUpdateTimeSteps[i] < InputTimeSteps[0] || 
           OutputUpdateTimeSteps[i] > InputTimeSteps[InputTimeStepsLength-1])
        {
            ERROR(<< "Requested Update Time Step should be in Time Range")
            vtkErrorMacro("Requested Update Time Step should be in Time Range.");
            return false;
        }

        // Loop over InputDataTimeSteps
        for(unsigned int InputDataTimeStepsIterator=0;
            InputDataTimeStepsIterator < InputDataTimeStepsLength;
            InputDataTimeStepsIterator++)
        {
            // Requested times that Snapped to data times do not need interpolation
            double SnapToTimeStepDifference = OutputUpdateTimeSteps[OutputUpdateTimeStepsIterator] - 
                InputTimeSteps[InputTimeStepsIterator];
            if(fabs(SnapToTimeiStepDifference) < this->SnapToTimeStepTolerance)
            {
                RequestedDataTimeIntervals[OutputUpdateTimeStepsIterator] = InputDataTimeStepsIterator;
                InterpolationRequested[OutputUpdateTimeStepsIterator] = false;
                IntervalFound = true;
                break;
            }

            // Find Interval
            else if(InputDataTimeStepsIterator < InputDataTimeStepsLength-1 &&
                    OutputUpdateTimeSteps[OutputUpdateTimeStepsIterator] > InputDataTimeSteps[InputDataTimeStepsIterator] && 
                    OutputUpdateTimeSteps[OutputUpdateTimeStepsIterator] < InputDataTimeSteps[InputDataTimeStepsIterator+1])
            {
                RequestedDataTimeIntervals[OutputUpdateTimeStepsiterator] = InputDataTimeStepsIterator;
                InterpolationRequested[OutputUpdateTimeStepsIterator] = true;
                IntervalFound = true;
                break;
            }
        }

        if(IntervalFound == false)
        {
            ERROR(<< "Interval not found.");
            vtkErrorMacro("Interval not found.");
            return false;
        }
    }

    DISPLAY(RequestedDataTimeIntervals,OutputUpdateTimeStepsLength);
    DISPLAY(InterpolationRequested,OutputUpdateTimeStepsLength);

    return true;
}

// ============================
// Convert Intervals To Indices
// ============================

bool Interpolator::ConvertIntervalsToIndices(
        unsigned int *Intervals,
        unsigned int NumberOfIntervals,
        unsigned int *Indices,            // Output
        unsigned int NumberOfIndices)     // Output
{
    // Indices Vector
    vtkstd::vector<unsigned int> IndicesVector;
    unsigned int IndicesVectorLength = 0;

    for(unsigned int IntervalsIterator=0; IntervalsIterator < NumberOfIntervals; IntervalsIterator++)
    {
        ExtendedRequestIndex.push_back(RequestIndex[i]);
        NumberOfExtendedIndices++;
        if(RequestInterpolation[i] == true)
        {
            ExtendedRequestIndex.push_back(RequestIndex[i] + 1);
            NumberOfExtendedIndices++;
        }
    }

//    DISPLAY(ExtendedRequestIndex,NumberOfExtendedIndices);

    // Remove repetitive indices
    std::vector<unsigned int> SqueezedRequestIndex;
    SqueezedRequestIndex.push_back(ExtendedRequestIndex[0]);
    unsigned int NumberOfSqueezedIndices = 1;

    if(NumberOfExtendedIndices > 1)
    {
        bool IndexRepeated = false;
        for(unsigned int i=1; i<NumberOfExtendedIndices; i++, IndexRepeated=false)
        {
            for(unsigned int j=0; j<i; j++)
            {
                if(ExtendedRequestIndex[j] == ExtendedRequestIndex[i])
                {
                    IndexRepeated = true;
                    break;
                }
            }
            if(IndexRepeated == false)
            {
//                ExtendedRequestIndex[NumberOfSqueezedIndices] = ExtendedRequestIndex[i];
                SqueezedRequestIndex.push_back(ExtendedRequestIndex[i]);
                NumberOfSqueezedIndices++;
            }
        }
    }

}

// ====
// Test
// ====

void Interpolator::Test()
{
    vtkExecutive *Executive = this->GetExecutive();
    vtkInformation *OutputInformation = Executive->GetOutputInformation(0);
    vtkDataObject *OutputDataObject = OutputInformation->Get(vtkDataObject::DATA_OBJECT());
    std::cout << "Interpolator: OutputDataObject: " << OutputDataObject << std::endl;
    if(OutputDataObject == NULL)
    {
        ERROR(<< "OutputDataObject is NULL");
    }
    else
    {
        vtkMultiBlockDataSet *OutputMultiBlockData = vtkMultiBlockDataSet::SafeDownCast(OutputDataObject);
        std::cout << "Interpolator: OutputMultiBlockData: " << OutputMultiBlockData << std::endl;
        if(OutputMultiBlockData == NULL)
        {
            ERROR(<< "OutputMultiBlockData is NULL.");
        }
        else
        {
            std::cout << "NEW: Number of Blocks: " << OutputMultiBlockData->GetNumberOfBlocks() << std::endl;
            vtkPolyData *OutputPolyData = vtkPolyData::SafeDownCast(OutputMultiBlockData->GetBlock(0));
            if(OutputPolyData == NULL)
            {
                ERROR(<< "OutputPolyData is NULL");
            }
            else
            {
                std::cout << "New: Number of Points: " << OutputPolyData->GetNumberOfPoints() << std::endl;
            }
        }
    }
}

// ============
// Request Data
// ============

int Interpolator::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet *input = vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Check input
    if(input == NULL)
    {
        ERROR(<< "input is NULL");
        vtkErrorMacro("input is NULL");
        return 0;
    }
    else if(input->GetNumberOfBlocks() < 1)
    {
        ERROR(<< "Number of blocks is zero.")
        vtkErrorMacro("Number of blocks is zero.")
        return 0;
    }
    else if(input->GetNumberOfPoints() < 1)
    {
        ERROR(<< "Number of input points is zero.");
        vtkErrorMacro("Number of input points is zero.");
        return 0;
    }

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    DISPLAY(output->GetNumberOfPoints());

    /// Test 1 ///
    vtkDataObject *OutputDataObject = this->GetOutputDataObject(0);
    vtkMultiBlockDataSet *OutputMultiBlock = vtkMultiBlockDataSet::SafeDownCast(OutputDataObject);
    if(!OutputDataObject)
    {
        ERROR(<< "OutputDataObject is NULL");
    }
    else
    {
        std::cout << "Interpoltor, OutputDataObject: " << OutputDataObject << std::endl;
    }
    if(!OutputMultiBlock)
    {
        ERROR(<< "OutpuMultiBlock is NULL");
    }
    else
    {
        std::cout << "Number of Blocks: " << OutputMultiBlock->GetNumberOfBlocks() << std::endl;
    }
    vtkPolyData *OutputPolyData = vtkPolyData::SafeDownCast(OutputMultiBlock->GetBlock(0));
    if(!OutputPolyData)
    {
        ERROR(<< "OutputPolyData is NULL.")
    }
    else
    {
        std::cout << "Number of points: " << OutputPolyData->GetNumberOfPoints() << std::endl;
    }

    // Check output
    if(output == NULL)
    {
        ERROR(<< "Output PolyData is NULL.");
        vtkErrorMacro("output PolyData is NULL");
        return 0;
    }
    else if(output->GetNumberOfPoints() < 1)
    {
        // Test
        WARNING(<< "Interpolator has no points.");
        vtkDataObject *OutputBlock = output->GetBlock(0);
        OutputBlock = NULL;

        // No point to Interpolate
        return 1;
    }

    // Time Steps
    unsigned int TimeStepsLength = inputInfo->Length(FilterInformation::UPDATE_TIME_STEPS());
    double *TimeSteps = inputInfo->Get(FilterInformation::TIME_STEPS());

    // Update Time Steps
    unsigned int UpdateTimeStepsLength = outputInfo->Length(FilterInformation::UPDATE_TIME_STEPS());
    double *UpdateTimeSteps = outputInfo->Get(FilterInformation::UPDATE_TIME_STEPS());

    // Interpolate
    bool InterpolateStatus = this->Interpolate(
            input,output,
            TimeSteps,TimeStepsLength,
            UpdateTimeSteps,UpdateTimeStepsLength);

    output->ShallowCopy(vtkDataSet::SafeDownCast(input->GetBlock(0)));

    HERE
    this->Test();

    return 1;
}

// ===========
// Interpolate
// ===========

bool Interpolator::Interpolate(
        vtkMultiBlockDataSet *InputMultiBlockData,
        vtkMultiBlockDataSet *OutputMultiBlockData,
        double *TimeSteps,
        unsigned int TimeStepsLength,
        double *UpdateTimeSteps,
        unsigned int UpdateTimeStepsLength)
{
    // 1- Inputs //

    // Number of Input Blocks or DataObjects
    unsigned int NumberOfInputDataObjects = InputMultiBlockData->GetNumberOfBlocks();


    // Using DataObject Array instead of MultiBlock Data
    vtkDataObject *InputDataObjectsArray[NumberOfInputDataObjects];

    // Convert MultiBlock Data to DataObject Array
    this->MultiBlockDataSetToDataObjectArray(InputMultiBlockData,InputDataObjectsArray);

    // 2- Outputs //
   
    // Number of Output Blocks or DataObjects
    unsigned int NumberOfOutputDataObjects = OutputMultiBlockData->GetNumberOfBlocks();

    // Using DataObject Array instead of MultoBlock Data
    vtkDataObject *OutputDataObjectsArray[NumberOfOutputDataObjects];

    // Convert MultiBlock Data to DataObject Array
    this->MultiBlockDataSetToDataObjectArray(OutputMultiBlockData,OutputDataObjectsArray);

    // Declare arrays needed for interpolation //

    // 1- CellIds of Tracers inside Data grid
    vtkIdTypeArray  *TracersCellIdsArray[NumberOfOutputDataObjects];

    // 2- Local coordinates of tracers in the found cell
    vtkDoubleArray *TracerParametricCoordinatesInCellArray[NumberOfOutputDataObjects];

    // 3- Weights of tracers in the cell found
    vtkDoubleArray *TracerWeightsInCellArray[NumberOfOutputDataObjects];

    // Define arrays needed for intrpolation //

    // Loop over output blocks
    for(unsigned int OutputDataObjectsIterator = 0;
        OutputDataObjectsIterator < NumberOfOutputDataObjects;
        OutputDataObjectsIterator++)
    {
        // TODO : delete the folowings on memory after the work is done.

        // 1- CellIds of Tracers inside Data Grid
        TracersCellIdsArray[OutputDataObjectsIterator] = vtkIdTypeArray::New();

        // 2- Local coordinates of tracers in the found cell
        TracerParametricCoordinatesInCellArray[OutputDataObjectsIterator] = vtkDoubleArray::New();

        // 3- Weights of tracers in the cell found
        TracerWeightsInCellArray[OutputDataObjectsIterator] = vtkDoubleArray::New();
    }

    // Get the first of DataObject for points for geometry
    vtkDataObject *FirstInputDataObject = InputDataObjectsArray[0];

    // Check First Input DataObject
    if(FirstInputDataObject == NULL)
    {
        ERROR(<< "FirstInputDataObject is NULL.");
        vtkErrorMacro("FirstInputDataObject is NULL.");
    }
    vtkDataSet *FirstInputDataSet = vtkDataSet::SafeDownCast(FirstInputDataObject);
    if(FirstInputDataSet->GetNumberOfPoints() < 1)
    {
        ERROR(<< "FirstInputDataSet has no point.");
        vtkErrorMacro("FirstInputDataSet has no points.");
    }

    // Locate tracer points inside data grid
    this->LocateTracersInDataGrid(
            FirstInputDataObject,
            OutputDataObjectsArray,
            NumberOfOutputDataObjects,
            TracersCellIdsArray,                       // Output
            TracerParametricCoordinatesInCellArray,    // Output
            TracerWeightsInCellArray);                 // Output

    // Loop over UpdateTimeSteps
    for(unsigned int UpdateTimeStepsIterator = 0;
        UpdateTimeStepsIterator < UpdateTimeStepsLength;
        UpdateTimeStepsIterator++)
    {
        // Find UpdateTimeSteps Index among intervals of TimeSteps array
        int UpdateTimeStepIntervalIndex = this->FindIntervalIndexInArray(
                UpdateTimeSteps[UpdateTimeStepsIterator],TimeSteps,TimeStepsLength);

        // Check if UpdateTimeStep is inside TimeRange
        if(UpdateTimeStepIntervalIndex < 0)
        {
            ERROR(<< "UpdateTimeStep is out of range.");
            vtkErrorMacro("UpdateTimeStep is out of range.");
        }

        // Find temporal intepolation coefficient
        double TemporalInterpolationCoefficient = this->FindTemporalInterpolationCoefficient(
                UpdateTimeSteps[UpdateTimeStepsIterator],
                TimeSteps[UpdateTimeStepIntervalIndex],
                TimeSteps[UpdateTimeStepIntervalIndex+1]);

        // Loop over Output DataObject Array

        for(unsigned int OutputDataObjectsIterator = 0;
            OutputDataObjectsIterator < NumberOfOutputDataObjects;
            OutputDataObjectsIterator++)
        {
            /// TEST ///
            DISPLAY(UpdateTimeStepIntervalIndex)
            DISPLAY(TimeSteps,TimeStepsLength)
            DISPLAY(UpdateTimeSteps,UpdateTimeStepsLength)

            // Spatial Interpolation
            this->SpatioTemporalInterpolation(
                    InputDataObjectsArray[UpdateTimeStepIntervalIndex],
                    InputDataObjectsArray[UpdateTimeStepIntervalIndex+1],
                    OutputDataObjectsArray[OutputDataObjectsIterator],
                    TracersCellIdsArray[OutputDataObjectsIterator],
                    TracerWeightsInCellArray[OutputDataObjectsIterator],
                    TemporalInterpolationCoefficient);
        }
    }

    // Append Tracer's CellIds to Output DataObjects
    for(unsigned int OutputDataObjectsIterator = 0;
        OutputDataObjectsIterator < NumberOfOutputDataObjects;
        OutputDataObjectsIterator++)
    {
        // Cast DataObject to DataSet
        vtkDataSet *AnOutputDataSet = vtkDataSet::SafeDownCast(OutputDataObjectsArray[OutputDataObjectsIterator]);

        // Set scalars for DataSet
        AnOutputDataSet->GetPointData()->SetScalars(TracersCellIdsArray[OutputDataObjectsIterator]);
    }

    // Freed memory
    for(unsigned int OutputDataObjectsIterator = 0;
        OutputDataObjectsIterator < NumberOfOutputDataObjects;
        OutputDataObjectsIterator++)
    {
        TracerParametricCoordinatesInCellArray[OutputDataObjectsIterator]->Delete();
        TracerWeightsInCellArray[OutputDataObjectsIterator]->Delete();
    }
    HERE

    return true;
}

// ===================================
// MultiBlock Data to DataObject Array
// ===================================

void Interpolator::MultiBlockDataSetToDataObjectArray(
        vtkMultiBlockDataSet *MultiBlockData,
        vtkDataObject **DataObjectArray)
{
    // Check MultiBlock
    if(MultiBlockData == NULL)
    {
        ERROR(<< "MultiBlock Data is NULL");
        vtkErrorMacro("MultiBlock Data is NULL.");
    }

    // Number of Blocks
    unsigned int NumberOfBlocks = MultiBlockData->GetNumberOfBlocks();

    if(NumberOfBlocks == 0)
    {
        ERROR(<< "MultiBlock has no leaf.");
        vtkErrorMacro("MultiBlock has no leaf.");
    }

    // Iterate over blocks
    for(unsigned int DataObjectIterator = 0; DataObjectIterator < NumberOfBlocks; DataObjectIterator++)
    {
        // Get DataObject block
        DataObjectArray[DataObjectIterator] = MultiBlockData->GetBlock(DataObjectIterator);

        // Check Data is valid
        if(DataObjectArray[DataObjectIterator] == NULL)
        {
            ERROR(<< "DataObject is NULL.");
            vtkErrorMacro("DataObject is NULL.");
        }
        else
        {
            // Cast to DataSet
            vtkDataSet *DataSet = vtkDataSet::SafeDownCast(DataObjectArray[DataObjectIterator]);

            // Check if number of points are zero
            if(DataSet->GetNumberOfPoints() < 1)
            {
                ERROR(<< "DataSet has zero number of points.")
                vtkErrorMacro("DataSet has zero number of points.");
            }
        }
    }
}

// ============================
// Find Interval Index In Array
// ============================

int Interpolator::FindIntervalIndexInArray(
        double InquiryValue,
        double *Array,
        unsigned int ArrayLength)
{
    /// TEST ///
    DISPLAY(InquiryValue);
    DISPLAY(Array,ArrayLength);
    // Initialize out of range flag
    int OUT_OF_RANGE_FLAG = -1;

    // Interval Index
    int IntervalIndex = OUT_OF_RANGE_FLAG;

    // Loop over Intervals of array
    for(unsigned int i=0; i<ArrayLength-1; i++)
    {
        // Check if value is in the interval
        if(InquiryValue >= Array[i] && InquiryValue <= Array[i+1])
        {
            IntervalIndex = i;
            break;
        }
    }

    // return -1 if no interval found
    return IntervalIndex;
}

// =======================================
// Find Temporal Interpolation Coefficient
// =======================================

double Interpolator::FindTemporalInterpolationCoefficient(
        double UpdateTimeStep,
        double TimeStepLeft,
        double TimeStepRight)
{
    // Check if time is not out of interval range
    if( UpdateTimeStep > TimeStepRight || UpdateTimeStep < TimeStepLeft)
    {
        ERROR(<< "UpdateTimeStep is out of interval range.");
        vtkErrorMacro("UpdateTimeStep is out of interval range.");
    }

    // Interval Range
    double IntervalRange = TimeStepRight - TimeStepLeft;

    // Interpolator Coefficient
    double InterpolationCoefficient = (UpdateTimeStep-TimeStepLeft) / IntervalRange;

    return InterpolationCoefficient;
}

// ===========================
// Locate Tracers in Data Grid
// ===========================

void Interpolator::LocateTracersInDataGrid(
        vtkDataObject *FirstInputDataObject,
        vtkDataObject **OutputDataObjectsArray,
        unsigned int NumberOfOutputDataObjects,
        vtkIdTypeArray **TracersCellIdsArray,
        vtkDoubleArray **TracerParametricCoordinatesInCellArray,
        vtkDoubleArray **TracerWeightsInCellArray)
{
    // Convert Input DataGrid to DataSet
    vtkDataSet *FirstInputDataSet = vtkDataSet::SafeDownCast(FirstInputDataObject);

    // Check FirstInputDataSet
    if(FirstInputDataSet == NULL)
    {
        ERROR(<< "First Input DataSet is NULL.");
        vtkErrorMacro("First Input DataSet is NULL.");
    }
    else if(FirstInputDataSet->GetNumberOfCells() < 1)
    {
        ERROR(<< "First Input DataSet does not have any cell.");
        vtkErrorMacro("First Input DataSet does not have any cell.");
    }

    // Get number of points of cells in Input DataSet
    unsigned int NumberOfPointsInACell = FirstInputDataSet->GetCell(1)->GetNumberOfPoints();

    // Loop over Output DataObjects (Blocks)
    for(unsigned int OutputDataObjectsIterator = 0;
        OutputDataObjectsIterator < NumberOfOutputDataObjects;
        OutputDataObjectsIterator++)
    {
        // Cast DataObject to DataSet
        vtkDataSet *OutputDataSet = vtkDataSet::SafeDownCast(OutputDataObjectsArray[OutputDataObjectsIterator]);

        // Number of Tracers
        unsigned int NumberOfTracers = OutputDataSet->GetNumberOfPoints();

        // Set Tracer CellIds
        TracersCellIdsArray[OutputDataObjectsIterator]->SetNumberOfComponents(1);
        TracersCellIdsArray[OutputDataObjectsIterator]->SetNumberOfTuples(NumberOfTracers);
        TracersCellIdsArray[OutputDataObjectsIterator]->SetName("CellIds");

        // Set Tracer Coordinates In Cell
        TracerParametricCoordinatesInCellArray[OutputDataObjectsIterator]->SetNumberOfComponents(this->Dimension);
        TracerParametricCoordinatesInCellArray[OutputDataObjectsIterator]->SetNumberOfTuples(NumberOfTracers);
        TracerParametricCoordinatesInCellArray[OutputDataObjectsIterator]->SetName("TracerParametricCoordinates");

        // Set Tracer Weights In Cell
        TracerWeightsInCellArray[OutputDataObjectsIterator]->SetNumberOfComponents(NumberOfPointsInACell);
        TracerWeightsInCellArray[OutputDataObjectsIterator]->SetNumberOfTuples(NumberOfTracers);
        TracerWeightsInCellArray[OutputDataObjectsIterator]->SetName("TracerWeights");

        // Get Tracers Previous cell Ids
        vtkIdTypeArray *TracersPreviousCellIds = vtkIdTypeArray::SafeDownCast(OutputDataSet->GetPointData()->GetScalars("CellIds"));

        // Check if Previous cell Ids exists
        bool TracersPreviousCellIdsExists = true;

        if(TracersPreviousCellIds == NULL)
        {
            TracersPreviousCellIdsExists = false;
        }
        else if(TracersPreviousCellIds->GetNumberOfTuples() != NumberOfTracers)
        {
            TracersPreviousCellIdsExists = false;
        }

        // Locate points in DataGrid //
        
        // Tracer point for inquiry
        double TracerPoint[this->Dimension];

        // Tracer Cell Id in Input Data grid
        vtkIdType TracerCellId;

        // Variables used for vtkDataSet::FindCellA
        vtkCell *GuessCell = NULL;
        vtkIdType GuessCellId = 0;
        double Tolerance = EPSILON;
        int TempSubId;

        // Output variables of vtkDataSet::FindCell
        double TracerParametricCoordinates[this->Dimension];
        double TracerWeights[NumberOfPointsInACell];

        // Iterate over all tracer points in current block
        for(unsigned int TracerIterator = 0;
            TracerIterator < NumberOfTracers;
            TracerIterator++)
        {
            // Get tracer point
            OutputDataSet->GetPoint(TracerIterator,TracerPoint);

            // Use local serach if a guess cell exists
            if(TracersPreviousCellIdsExists)
            {
                // Set guess cell from previous cell id
                GuessCellId = TracersPreviousCellIds->GetValue(TracerIterator);
            }

            // Find cell
            if(GuessCellId >= 0)
            {
                // Previous point was in domain
                TracerCellId = FirstInputDataSet->FindCell(
                        TracerPoint,
                        GuessCell,
                        GuessCellId,
                        Tolerance,
                        TempSubId,
                        TracerParametricCoordinates,   // Output
                        TracerWeights);                // Output
            }
            else
            {
                // Previous point was not in domain
                TracerCellId = -1;
            }

            // Set CellId in array
            TracersCellIdsArray[OutputDataObjectsIterator]->SetValue(TracerIterator,TracerCellId);

            // Set Tracer coordinates in array
            TracerParametricCoordinatesInCellArray[OutputDataObjectsIterator]->SetTuple(TracerIterator,TracerParametricCoordinates);

            // Set tracer weights in array
            TracerWeightsInCellArray[OutputDataObjectsIterator]->SetTuple(TracerIterator,TracerWeights);
        }        
    }
}

// =============================
// Spatio-Temporal Interpolation
// =============================

void Interpolator::SpatioTemporalInterpolation(
        vtkDataObject *InputDataObjectLeft,
        vtkDataObject *InputDataObjectRight,
        vtkDataObject *OutputDataObject,
        vtkIdTypeArray *TracerCellIds,
        vtkDoubleArray *TracerWeightsInCell,
        double TemporalInterpolationCoefficient)
{
    // 1- On the left of time interval //

    // Input DataSet on Left of time intrval
    vtkDataSet *InputDataSetLeft = vtkDataSet::SafeDownCast(InputDataObjectLeft);

    // Check Input DataSet Left
    if(InputDataSetLeft == NULL)
    {
        ERROR(<< "InputDataSetLeft is NULL.");
        vtkErrorMacro("InputDataSetLeft is NULL.");
    }
    else if(InputDataSetLeft->GetNumberOfPoints() < 1)
    {
        ERROR(<< "InputDataSetLeft has no points.");
        vtkErrorMacro("InputDataSetLeft has no points.");
    }

    // Velocities on Left of time interval
    vtkDoubleArray *InputVelocitiesLeft = vtkDoubleArray::SafeDownCast(
            InputDataSetLeft->GetPointData()->GetVectors("velocity"));
            // InputDataSetLeft->GetPointData()->GetVectors("DataVelocities"));

    // Check Left Input Velocities
    if(InputVelocitiesLeft == NULL)
    {
        ERROR(<< "Input velocities for left of time interval is NULL.");
        vtkErrorMacro("Input velocities for left of time interval is NULL");
    }
    else if(InputVelocitiesLeft->GetNumberOfTuples() < 1)
    {
        ERROR(<< "Velocities of input data on the left of time interval has no tuples.");
        vtkErrorMacro("Velocities of input data on the left of time interval has no tuples.");
    }
    else if(InputVelocitiesLeft->GetNumberOfComponents() != this->Dimension)
    {
        ERROR(<< "Dimension of velocity data are not consistent.");
        vtkErrorMacro("Dimension of velocity data are not consistent.");
    }
    HERE
    InputDataObjectRight->Print(std::cout);

    // 2- On the Right of time interval //

    // Input DataSet on Right of time interval
    vtkDataSet *InputDataSetRight = vtkDataSet::SafeDownCast(InputDataObjectRight);
    HERE

    // Check Input DataSet Right
    if(InputDataSetRight == NULL)
    {
        ERROR(<< "InputDataSetRight is NULL.");
        vtkErrorMacro("InputDataSetRight is NULL.");
    }
    else if(InputDataSetRight->GetNumberOfPoints() < 1)
    {
        ERROR(<< "InputDataSetRight has no points.");
        vtkErrorMacro("InputDataSetRight has no points.");
    }
    HERE

    // Velocities data on Right of time interval
    vtkDoubleArray *InputVelocitiesRight = vtkDoubleArray::SafeDownCast(
            InputDataSetRight->GetPointData()->GetVectors("velocity"));
    HERE

    // Check Right Input Velocities
    if(InputVelocitiesLeft == NULL)
    {
        ERROR(<< "Input velocities for left of time interval is NULL.");
        vtkErrorMacro("Input velocities for left of time interval is NULL");
    }
    else if(InputVelocitiesLeft->GetNumberOfTuples() < 1)
    {
        ERROR(<< "Velocities of input data on the left of time interval has no tuples.");
        vtkErrorMacro("Velocities of input data on the left of time interval has no tuples.");
    }
    else if(InputVelocitiesLeft->GetNumberOfComponents() != this->Dimension)
    {
        ERROR(<< "Dimension of velocity data are not consistent.");
        vtkErrorMacro("Dimension of velocity data are not consistent.");
    }
    
    // 3- Output data //

    // Output DataSet
    vtkDataSet *OutputDataSet = vtkDataSet::SafeDownCast(OutputDataObject);

    // Number of Tracers
    unsigned int NumberOfTracers = OutputDataSet->GetNumberOfPoints();

    // Output (Tracers) velocities
    vtkDoubleArray *OutputVelocities = vtkDoubleArray::SafeDownCast(
            OutputDataSet->GetPointData()->GetVectors("TracerVelocities"));

    // Define new velocity array
    if(OutputVelocities == NULL)
    {
        OutputVelocities = vtkDoubleArray::New();
    }

    // Set Number of Tuples
    if(OutputVelocities->GetNumberOfTuples() != NumberOfTracers)
    {
        OutputVelocities->SetNumberOfTuples(NumberOfTracers);
    }

    // Set Number of Components
    if(OutputVelocities->GetNumberOfComponents() != this->Dimension)
    {
        OutputVelocities->SetNumberOfComponents(this->Dimension);
    }

    // Set array name again
    OutputVelocities->SetName("TracerVelocities");

    // Temporal Interpolatoin Coefficients Left and Right
    double TemporalInterpolationCoefficientLeft = 1 - TemporalInterpolationCoefficient;
    double TemporalInterpolationCoefficientRight = TemporalInterpolationCoefficient;

    // Interpolate for each tracer point
    for(unsigned int TracerIterator = 0; TracerIterator < NumberOfTracers; TracerIterator++)
    {
        // Get Tracer's CellId
        double TracerCellIdDoubleType = TracerCellIds->GetTuple1(TracerIterator);
        int TracerCellId = static_cast<int>(round(TracerCellIdDoubleType));

        // skip points outside data grid
        if(TracerCellId < 0)
        {
            // point is outside data grid
            continue;
        }

        // Get Tracer's Cell (Either from left or right data)
        vtkCell *TracerCell = InputDataSetLeft->GetCell(TracerCellId);

        // Get Number of points in the cell
        unsigned int NumberOfPointsInCell = TracerCell->GetNumberOfPoints();

        // Get Points in Cell
        vtkIdList *CellPointIds = TracerCell->GetPointIds();

        // Get Cell points weights
        double *CellPointsWeights = TracerWeightsInCell->GetTuple(TracerIterator);

        // Declare Tracer and cell point velocity
        double TracerPointVelocity[this->Dimension];
        double CellPointVelocityLeft[this->Dimension];
        double CellPointVelocityRight[this->Dimension];

        // Initialize tracer point velocity to zero
        for(unsigned int DimensionIterator = 0;
            DimensionIterator < this->Dimension;
            DimensionIterator++)
        {
            TracerPointVelocity[DimensionIterator] = 0;
        }

        // Loop over points of cell
        for(unsigned int CellPointsIterator = 0; CellPointsIterator < NumberOfPointsInCell; CellPointsIterator++)
        {
            // Get Cell Point Id
            vtkIdType CellPointId = CellPointIds->GetId(CellPointsIterator);

            // Get Cell points velocity
            InputVelocitiesLeft->GetTuple(CellPointId,CellPointVelocityLeft);
            InputVelocitiesRight->GetTuple(CellPointId,CellPointVelocityRight);

            // Calculate Output velocities, Iterate over dimension
            for(unsigned int DimensionIterator = 0;
                DimensionIterator < this->Dimension;
                DimensionIterator++)
            {
                TracerPointVelocity[DimensionIterator] += 
                    CellPointVelocityLeft[DimensionIterator] * CellPointsWeights[CellPointsIterator] * TemporalInterpolationCoefficientLeft +
                    CellPointVelocityRight[DimensionIterator] * CellPointsWeights[CellPointsIterator] * TemporalInterpolationCoefficientRight;
            }
        }
    }
}

// =====================
// Find Grid Adjacencies
// =====================

void Interpolator::FindGridAdjacencies(
        vtkDataObject *DataObjectGrid,
        vtkIdTypeArray *GridAdjacenciesArray)   // Output
{
    // Cast data object to data set
    vtkDataSet *DataSetGrid = vtkDataSet::SafeDownCast(DataObjectGrid);

    // Check Data set
    if(DataSetGrid == NULL)
    {
        ERROR(<< "DataSet Gird is NULL.");
        vtkErrorMacro("DataSet Grid is NULL.");
    }
    
    // Get Number of Cells
    vtkIdType NumberOfCells = DataSetGrid->GetNumberOfCells();

    // Check number of cells
    if(NumberOfCells < 1)
    {
        ERROR(<< "Number of Cells is zero.");
        vtkErrorMacro("Number of Cells is zero.");
    }

    // Determine grid topology //

    // Initialize maximum number of faces among all cells
    unsigned int MaximumNumberOfCellFaces = 0;

    // Initialize Number of CellFaces Array
    vtkIdTypeArray *NumberOfCellFacesArray = NULL;

    // Fine grid topology
    GridType DataSetGridType = FindGridTopology(
            DataSetGrid,
            MaximumNumberOfCellFaces,     // Output
            NumberOfCellFacesArray);      // Output

    // Define Grid Adjacencies Array
    GridAdjacenciesArray = vtkIdTypeArray::New();
    GridAdjacenciesArray->SetNumberOfTuples(NumberOfCells);
    GridAdjacenciesArray->SetNumberOfValues(MaximumNumberOfCellFaces);

    // Find Adjacencies //

    // Iterate over cells
    for(vtkIdType CellIterator = 0; CellIterator < NumberOfCells; CellIterator++)
    {
        // Get a cell
        vtkGenericCell *GenericCell = vtkSmartPointer<vtkGenericCell>::New();
        DataSetGrid->GetCell(CellIterator,GenericCell);

        // Get number of faces of the cell
        unsigned int NumberOfCellFaces;

        if(DataSetGridType == UNIFORM_GRID)
        {
            // For Uniform Grids
            NumberOfCellFaces = MaximumNumberOfCellFaces;
        }
        else
        {
            // For Non-Uniform (Hybrid) Grids
            NumberOfCellFaces = NumberOfCellFacesArray->GetValue(CellIterator);
        }

        // Loop over faces of the cell
        for(unsigned int CellFaceIterator = 0; CellFaceIterator < NumberOfCellFaces; CellFaceIterator++)
        {
            // Get a face
            vtkCell *CellFaceCell = GenericCell->GetFace(CellFaceIterator);

            // Get point Ids of cell that constitute the face
            vtkIdList *CellFaceCellPointIds = CellFaceCell->GetPointIds();

            // Define id lists of neighbor cells
            vtkIdList *NeighborCellIds = vtkSmartPointer<vtkIdList>::New();

            // Look for another cell sharing the same face
            DataSetGrid->GetCellNeighbors(CellIterator,CellFaceCellPointIds,NeighborCellIds);

            // Number of neighbor cells
            unsigned int NumberOfNeighborCells = NeighborCellIds->GetNumberOfIds();

            if(NumberOfNeighborCells > 0)
            {
                // Multiple neighbors
                if(NumberOfNeighborCells > 1)
                {
                    WARNING(<< "Ambigous number of cell neighbors.");
                }

                // Get neighbor cell Id
                vtkIdType NeighborCellId = NeighborCellIds->GetId(0);

                // Set adjacencies
                GridAdjacenciesArray->SetComponent(CellIterator,CellFaceIterator,NeighborCellId);
            }
            else
            {
                // No neighbor found, face is on grid boundary
                GridAdjacenciesArray->SetComponent(CellIterator,CellFaceIterator,-1);
            }
        }
    }
}

// ==================
// Find Grid Topology
// ==================

GridType Interpolator::FindGridTopology(
        vtkDataObject *DataObjectGrid,
        unsigned int &MaximumNumberOfCellFaces,      // Output
        vtkIdTypeArray *NumberOfCellFacesArray)      // Output
{
    // Cast to DataSet
    vtkDataSet *DataSetGrid = vtkDataSet::SafeDownCast(DataObjectGrid);

    // Check DataSet
    if(DataSetGrid == NULL)
    {
        ERROR(<< "DataSet is NULL.");
        vtkErrorMacro("DataSet is NULL");
        return UNDEFINED_GRID;
    }

    // Number of Cells
    unsigned int NumberOfCells = DataSetGrid->GetNumberOfCells();

    // Check number of cells
    if(NumberOfCells < 1)
    {
        ERROR(<< "Number of cells is zero.");
        vtkErrorMacro("NumberOfCells is zero.");
        return UNDEFINED_GRID;
    }

    // Define output array
    NumberOfCellFacesArray = vtkIdTypeArray::New();
    NumberOfCellFacesArray->SetNumberOfTuples(NumberOfCells);
    NumberOfCellFacesArray->SetNumberOfValues(1);

    // Assume grid type is unifirm unless otherwise detected
    GridType DataGridType = UNIFORM_GRID;

    // Initialize Maximum number of faces in cell from a sample cell
    MaximumNumberOfCellFaces = DataSetGrid->GetCell(0)->GetNumberOfFaces();

    // Loop over cells
    for(unsigned int CellIterator = 0; CellIterator < NumberOfCells; CellIterator++)
    {
        // Get a cell
        vtkGenericCell *GenericCell = vtkSmartPointer<vtkGenericCell>::New();
        DataSetGrid->GetCell(CellIterator,GenericCell);

        // Number of faces of cell
        unsigned int NumberOfCellFaces = GenericCell->GetNumberOfFaces();

        // Compare to maximum number of faces
        if(NumberOfCellFaces != MaximumNumberOfCellFaces)
        {
            // Change type of grid
            DataGridType = HYBRID_GRID;

            // Compare maximum number of faces
            if(NumberOfCellFaces > MaximumNumberOfCellFaces)
            {
                // Update maximum number of faces
                MaximumNumberOfCellFaces = NumberOfCellFaces;
            }
        }

        // Set Number of CellFaces array
        NumberOfCellFacesArray->SetValue(CellIterator,NumberOfCellFaces);
    }

    // If the grid is uniform, array is not needed
    if(DataGridType == UNIFORM_GRID)
    {
        // Delete array
        NumberOfCellFacesArray->Delete();
        NumberOfCellFacesArray == NULL;
    }

    // Output grid type
    return DataGridType;
}

// ===========================
// Find Cell with Local Search
// ===========================

vtkIdType Interpolator::FindCellWithLocalSearch(
        double *InquiryPoint,
        unsigned int Dimension,
        vtkIdType GuessCellId,
        vtkDataObject *DataObjectGrid,
        double *FoundCellWeights,            // Output
        unsigned int NumberOfWeights)        // Output
{
    // Cast DataObject to DataSet
    vtkDataSet *DataSetGrid = vtkDataSet::SafeDownCast(DataObjectGrid);

    // Get Guess Cell
    vtkGenericCell *GuessGenericCell = vtkSmartPointer<vtkGenericCell>::New();
    DataSetGrid->GetCell(GuessCellId,GuessGenericCell);

    // Get Cell Point Ids
    vtkIdList *CellPointIds = GuessGenericCell->GetPointIds();

    

}

// ===========================
// Triangle Cell Interpolation
// ===========================

bool Interpolator::TriangleCellInterpolation(
        const double *InquiryPoint,
        const unsigned int Dimension,
        const double **CellPoints,
        double *InquiryPointLocalCoordinates,    // Output
        unsigned int NumberOfLocalCoordinates,   // Output
        double *CellPointWights,                 // Output
        unsigned int NumberOfWeights)            // Output
{

}

// ==============================
// Tetrahedron Cell Interpolation
// ==============================

bool Interpolator::TetrahedronCellInterpolation(
        const double *InquiryPoint,
        const unsigned int Dimension,
        vtkPoints *CellPoints,
        double *InquiryPointLocallCoordinates,   // Output
        unsigned int NumberOfLocalCoordinates,   // Output
        double *CellPointWeights,                // Output
        unsigned int NumberOfWeights)            // Output
{
    // Number of coordinates
    NumberOfLocalCoordinates = 3;

    // Base point for coordinate
    double *BasePoint = CellPoints->GetPoint(0);

    // Edge vectors
    double *EdgeVectors[Dimension];

    // Set edge vectors
    for(unsigned int EdgeIterator = 0; EdgeIterator < Dimension; EdgeIterator++)
    {
        // Get Edge point
        double *EdgePoint = CellPoints->GetPoint(EdgeIterator+1);

        // Construct edge vector
        for(unsigned int DimensionIterator = 0; DimensionIterator < Dimension; DimensionIterator++)
        {
            EdgeVectors[EdgeIterator][DimensionIterator] = EdgePoint[DimensionIterator] - BasePoint[DimensionIterator];
        }
    }


}
