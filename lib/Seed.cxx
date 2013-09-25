/*
 * =====================================================================================
 *
 *       Filename:  Seed.cxx
 *
 *    Description:  Search for appropriate Seed points
 *
 *        Version:  1.0
 *        Created:  10/16/2012 09:12:57 AM
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

#include <Seed.h>

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

// Data
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>

// Elements
#include <vtkPoints.h>
#include <vtkCell.h>

// Arrays
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

// Keys
#include <vtkInformationIntegerKey.h>

// For DEBUG
#include <vtkInformationRequestKey.h>

// ======
// Macros
// ======

#define EPSILON 1e-4  // for search cells
vtkStandardNewMacro(Seed);
vtkCxxRevisionMacro(Seed,"$Revision 1.0$");

// ===========
// Constructor
// ===========

Seed::Seed()
{
    // Member Data
    this->Resolution = NULL;
    this->Bounds = NULL;
    this->Spacing = NULL;
    this->Origin = NULL;
    this->Dimension = 3;
    this->SearchMode = SEARCH_MODE_USE_vtkPointSet;

    // Pipeline
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

// ==========
// Destructor
// ==========

Seed::~Seed()
{
    delete [] this->Resolution;
    this->Resolution = 0;

    delete [] this->Bounds;
    this->Bounds = 0;

    delete [] this->Spacing;
    this->Spacing = 0;

    delete [] this->Origin;
    this->Origin = 0;
}

// ==========
// Print Self
// ==========

void Seed::PrintSelf(ostream &os,vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ===================
// Accessors, Mutators
// ===================

// Resolution
int * Seed::GetResolution() const
{
    return this->Resolution;
}

void Seed::GetResolution(int *resolution)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.");
        vtkErrorMacro("Dimension should be defined first.");
    }

    resolution = this->Resolution;
}

void Seed::SetResolution(int *resolution)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.");
        vtkErrorMacro("Dimension should be defined first.");
    }

    this->Resolution = new int[this->Dimension];
    for(unsigned int i=0; i<this->Dimension; i++)
    {
        this->Resolution[i] = resolution[i];
    }
}

void Seed::SetResolution(int resolution0, int resolution1)
{
    if(this->Dimension != 2)
    {
        ERROR(<< "Dimension should be 2.");
        vtkErrorMacro("Dimension should be 2.");
    }

    this->Resolution = new int[this->Dimension];
    this->Resolution[0] = resolution0;
    this->Resolution[1] = resolution1;
}

void Seed::SetResolution(int resolution0, int resolution1, int resolution2)
{
    if(this->Dimension != 3)
    {
        ERROR(<< "Dimension should be 3.");
        vtkErrorMacro("Dimension should be 3.");
    }

    this->Resolution = new int[this->Dimension];
    this->Resolution[0] = resolution0;
    this->Resolution[1] = resolution1;
    this->Resolution[2] = resolution2;
}

// Bounds
double * Seed::GetBounds() const
{
    return this->Bounds;
}

void Seed::GetBounds(double *bounds)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.");
        vtkErrorMacro("Dimension should be defined first.");
    }

    bounds = this->Bounds;
}

void Seed::SetBounds(double *bounds)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.");
        vtkErrorMacro("Dimension should be defined first.");
    }

    this->Bounds = new double[2*this->Dimension];
    for(unsigned int i=0; i<2*this->Dimension; i++)
    {
        this->Bounds[i] = bounds[i];
    }
}

void Seed::SetBounds(double bounds0, double bounds1, double bounds2, double bounds3)
{
    if(this->Dimension !=2 )
    {
        ERROR(<< "Dimension should be 2.");
        vtkErrorMacro("Dimension should be 2.")
    }
    this->Bounds = new double[2*this->Dimension];
    this->Bounds[0] = bounds0;
    this->Bounds[1] = bounds1;
    this->Bounds[2] = bounds2;
    this->Bounds[3] = bounds3;
}

void Seed::SetBounds(double bounds0, double bounds1, double bounds2, double bounds3, double bounds4, double bounds5)
{
    if(this->Dimension !=3 )
    {
        ERROR(<< "Dimension should be 3.");
        vtkErrorMacro("Dimension should be 3.")
    }
    this->Bounds = new double[2*this->Dimension];
    this->Bounds[0] = bounds0;
    this->Bounds[1] = bounds1;
    this->Bounds[2] = bounds2;
    this->Bounds[3] = bounds3;
    this->Bounds[4] = bounds4;
    this->Bounds[5] = bounds5;
}

// Spacing
double * Seed::GetSpacing() const
{
    return this->Spacing;
}

void Seed::GetSpacing(double *spacing)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.")
            vtkErrorMacro("Dimension should be defined first.")
    }

    spacing = this->Spacing;
}

void Seed::SetSpacing(double *spacing)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.");
        vtkErrorMacro("Dimension should be defined first.")
    }
    this->Spacing = new double[this->Dimension];
    for(unsigned int i=0; i<this->Dimension; i++)
    {
        this->Spacing[i] = spacing[i];
    }
}

void Seed::SetSpacing(double spacing0, double spacing1)
{
    if(this->Dimension != 2)
    {
        ERROR(<< "Dimension should be 2.");
        vtkErrorMacro("Dimension should be 2.");
    }
    this->Spacing = new double[this->Dimension];
    this->Spacing[0] = spacing0;
    this->Spacing[1] = spacing1;
}

void Seed::SetSpacing(double spacing0, double spacing1,double spacing2)
{
    if(this->Dimension != 3)
    {
        ERROR(<< "Dimension should be 3.");
        vtkErrorMacro("Dimension should be 3.");
    }
    this->Spacing = new double[this->Dimension];
    this->Spacing[0] = spacing0;
    this->Spacing[1] = spacing1;
    this->Spacing[2] = spacing2;
}

// Origin
double * Seed::GetOrigin() const
{
    return this->Origin;
}

void Seed::GetOrigin(double *origin)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.");
        vtkErrorMacro("Dimension should be defined first.");
    }

    origin = this->Origin;
}

void Seed::SetOrigin(double *origin)
{
    if(this->Dimension < 2)
    {
        ERROR(<< "Dimension should be defined first.");
        vtkErrorMacro("Dimension should be defined first.");
    }

    this->Origin = new double[this->Dimension];
    for(unsigned int i=0; i<this->Dimension; i++)
    {
        this->Origin[i] = origin[i];
    }
}

void Seed::SetOrigin(double origin0, double origin1)
{
    if(this->Dimension != 2)
    {
        ERROR(<< "Dimension should be 2.");
        vtkErrorMacro("Dimension should be 2.");
    }

    this->Origin = new double[this->Dimension];
    this->Origin[0] = origin0;
    this->Origin[1] = origin1;
}

void Seed::SetOrigin(double origin0, double origin1, double origin2)
{
    if(this->Dimension != 3)
    {
        ERROR(<< "Dimension should be 3.");
        vtkErrorMacro("Dimension should be 3.");
    }

    this->Origin = new double[this->Dimension];
    this->Origin[0] = origin0;
    this->Origin[1] = origin1;
    this->Origin[2] = origin2;
}

// ==========
// Get Output
// ==========

vtkStructuredPoints * Seed::GetOutput()
{
    return this->GetOutput(0);
}

vtkStructuredPoints * Seed::GetOutput(int port)
{
    return vtkStructuredPoints::SafeDownCast(this->GetOutputDataObject(port));
}

// ==========
// Set Output
// ==========

void Seed::SetOutput(vtkDataObject *OutputDataObject)
{
    this->SetOutput(0,OutputDataObject);
}

void Seed::SetOutput(int port,vtkDataObject *OutputDataObject)
{
    this->GetExecutive()->SetOutputData(port, OutputDataObject);
}

// =========
// Get Input
// =========

vtkDataSet * Seed::GetInput()
{
    return this->GetInput(0);
}

vtkDataSet * Seed::GetInput(int port)
{
    return vtkDataSet::SafeDownCast(this->GetExecutive()->GetInputData(port,0));
}

// =========
// Set Input
// =========

void Seed::SetInput(vtkDataObject *InputDataObject)
{
    this->SetInput(0,InputDataObject);
}

void Seed::SetInput(int port, vtkDataObject *InputDataObject)
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

void Seed::AddInput(vtkDataObject *InputDataObject)
{
    this->AddInput(0,InputDataObject);
}

void Seed::AddInput(int port, vtkDataObject *InputDataObject)
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

int Seed::ProcessRequest(
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

int Seed::FillInputPortInformation(int port,vtkInformation *info)
{
    if(port == 0)
    {
        // TODO: set this for two different types with APPEND.
        // also change it in Seed.xml
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkDataObject");

        DEBUG(<< "Success");
        return 1;
    }

    DEBUG(<< "Failure");
    return 0;
}

// ============================
// Fill Output Port Information
// ============================

int Seed::FillOutputPortInformation(int port,vtkInformation *info)
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

int Seed::RequestDataObject(
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

int Seed::RequestInformation(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *vtkNotUsed(outputVector))
{
    // Input Info
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Check inputInfo is not NULL
    if(inputInfo == NULL)
    {
        ERROR(<< "inputInfo is NULL.");
        vtkErrorMacro("inputInfo is NULL.");
        return 0;
    }

    // 1- Time Steps //

    // Check if TimeSteps key exists in the input
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
        ERROR(<< "TimeSteps length is zero.");
        vtkErrorMacro("TimeSteps length is zero.");
        return 0;
    }

    if(TimeSteps == NULL)
    {
        ERROR(<< "TimeSteps is NULL.");
        vtkErrorMacro(<< "TimeSteps is NULL.");
        return 0;
    }

    // 2- Update Time Steps //

    // Create Update Time Steps
    unsigned int UpdateTimeStepsLength = 1;
    double UpdateTimeSteps[UpdateTimeStepsLength];
    
    // Set Update TimeSteps to the first time step
    UpdateTimeSteps[0] = TimeSteps[0];

    // Set UpdateTimeSteps key to input
    inputInfo->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);

    DEBUG(<< "Success");
    return 1;
}

// =====================
// Request Update Extent
// =====================

int Seed::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // input Info
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);

    // Check inputInfo is not NULL
    if(inputInfo == NULL)
    {
        ERROR(<< "inputInfo is NULL.");
        vtkErrorMacro("inputInfo is NULL.");
        return 0;
    }

    // 1- TIME STEPS //

    // Check inputInfo has TIME_STEPS key
    if(!inputInfo->Has(FilterInformation::TIME_STEPS()))
    {
        ERROR(<< "inputInfo does not have TIME_STEPS key.");
        vtkErrorMacro("inputInfo does not have TIME_STEPS key.");
        return 0;
    }

    // Get TimeSteps
    double *TimeSteps = inputInfo->Get(FilterInformation::TIME_STEPS());
    unsigned int TimeStepsLength = inputInfo->Length(FilterInformation::TIME_STEPS());

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

    // 2- UPDATE TIME STEPS //

    // Check inputInfo has UPDATE TIME STEPS key
    if(!inputInfo->Has(FilterInformation::UPDATE_TIME_STEPS()))
    {
        ERROR(<< "inputInfo doea not have UPDATE_TIME_STEPS key.");
        vtkErrorMacro("inputInfo does not have UPDATE_TIME_STEPS key.");
        return 0;
    }

    // Update Time Steps
    double *UpdateTimeSteps = &TimeSteps[0];
    unsigned int UpdateTimeStepsLength = 1;

    // Set UpdateTimeSteps
    inputInfo->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);

    DEBUG(<< "Success");
    return 1;
}

// ============
// Request Data
// ============

int Seed::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkDataObject *InputDataObject = inputInfo->Get(vtkDataObject::DATA_OBJECT());
    vtkDataSet *input = NULL;

    // Check inputInfo
    if(inputInfo == NULL)
    {
        ERROR(<< "inputInfo is NULL.");
        vtkErrorMacro("inputInfo is NULL.");
        return 0;
    }

    // Check Input DataObject
    if(InputDataObject == NULL)
    {
        ERROR(<< "InputDataObject is NULL.");
        vtkErrorMacro("InputDataObject is NULL.");
        return 0;
    }

    // Get one sample of MultiBlock

    // 1- If InputDataObject is MultiBlock
    if(strcmp(InputDataObject->GetClassName(),"vtkMultiBlockDataSet") == 0)
    {
        vtkMultiBlockDataSet *InputMultiBlockDataSet = vtkMultiBlockDataSet::SafeDownCast(InputDataObject);
        input = vtkDataSet::SafeDownCast(InputMultiBlockDataSet->GetBlock(0));
    }

    // 2- If InputDataObject is DataSet
    else
    {
        input = vtkDataSet::SafeDownCast(InputDataObject);
    }

    // Check Input DataSet
    if(input == NULL)
    {
        ERROR(<< "Input DataSet is NULL");
        vtkErrorMacro("Input DataSet is NULL.");
        return 0;
    }

    if(input->GetNumberOfPoints() <= 1)
    {
        vtkErrorMacro("Input data is empty or do not have enough points.");
        return 1;
    }

    // Output (StructuredGrid)
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkStructuredGrid *output = vtkStructuredGrid::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Output (StructuredPoints)
    vtkSmartPointer<vtkStructuredPoints> StructuredPointsOutput = vtkSmartPointer<vtkStructuredPoints>::New();

    // Set geometry of output grid
    this->SetOutputGridGeometry(StructuredPointsOutput);

    // Define CellIds array
    vtkSmartPointer<vtkIntArray> CellIds = vtkSmartPointer<vtkIntArray>::New();
    CellIds->SetName("CellIds");

    // Intersection of Grids
    FindIntersectionOfGrids(input,StructuredPointsOutput,CellIds);

    // Initialize Tracers
    InitializeTracers(input,StructuredPointsOutput,CellIds);

    // Add Array to Output
    output->GetPointData()->SetScalars(CellIds);

    // Convert StructuredPoints output to StructuredGrid output
    Seed::ConvertStructuredPointsToStructuredGrid(StructuredPointsOutput,output);

    DISPLAY(output->GetNumberOfPoints());
    DISPLAY(output->GetClassName());

    if(output->GetNumberOfPoints() > 0)
    {
        DEBUG(<< "Success");
        return 1;
    }
    else
    {
        DEBUG(<< "Failue");
        vtkErrorMacro("Number of points in output is zero.");
        return 0;
    }
}

// ========================
// Set Output Grid Geometry
// ========================

int Seed::SetOutputGridGeometry(vtkStructuredPoints *StructuredPoints)
{
    // Check geometry
    // 1- by Resolution and Bounds
    if(this->Resolution != NULL && this->Bounds != NULL)
    {
        this->Origin = new double[this->Dimension];
        this->Spacing = new double[this->Dimension];
        for(unsigned int i=0; i<this->Dimension; i++)
        {
            this->Origin[i] = this->Bounds[2*i];
            this->Spacing[i] = fabs(this->Bounds[2*i+1]-this->Bounds[2*i])/(this->Resolution[i]-1);
        }
    }

    // 2- by Bounds and Spacing
    else if(this->Bounds != NULL && this->Spacing != NULL)
    {
        this->Origin = new double[this->Dimension];
        this->Resolution = new int[this->Dimension];
        for(unsigned int i=0; i<this->Dimension; i++)
        {
            this->Origin[i] = this->Bounds[2*i];
            this->Resolution[i] = floor(fabs(this->Bounds[2*i+1]-this->Bounds[2*i])/this->Spacing[i]);
            if((fabs(this->Bounds[2*i+1]-this->Bounds[2*i])/this->Spacing[i]) - this->Resolution[i] > EPSILON)
            {
                WARNING(<< "Seed grid geometry is not exactly the same is requested extent.");
            }
        }
    }

    // 3- by Origin, Spacing and Resolution
    else if(this->Origin != NULL && this->Spacing != NULL && this->Resolution != NULL)
    {
        this->Bounds = new double[2*this->Dimension];
        for(unsigned int i=0; i<this->Dimension; i++)
        {
            this->Bounds[2*i] = this->Origin[i];
            this->Bounds[2*i+1] = this->Origin[i] + this->Spacing[i]*(this->Resolution[i]-1);
        }
    }

    // 4- Else, geometry undefined
    else
    {
        ERROR(<< "Not enough input data to define seed grid geometry.");
        vtkErrorMacro("Not enough input data for Seed grid. Geometry undefined.");
        return 0;
    }

    // Define Output Geometry
    if(this->Dimension == 3)
    {
        StructuredPoints->SetDimensions(this->GetResolution());
        StructuredPoints->SetOrigin(this->Origin);
        StructuredPoints->SetSpacing(this->Spacing);
    }
    else if(this->Dimension == 2)
    {
        StructuredPoints->SetDimensions(Resolution[0],Resolution[1],1);
        StructuredPoints->SetOrigin(this->Origin[0],this->Origin[1],0);
        StructuredPoints->SetSpacing(this->Spacing[0],this->Spacing[1],0);
    }
    else
    {
        vtkErrorMacro("Only 2 and 3 dimensions are supported.");
        return 0;
    }

    return 1;
}

// ==========================
// Find Intersection Of Grids
// ==========================

int Seed::FindIntersectionOfGrids(
        vtkDataSet *DataGrid,
        vtkDataSet *SeedGrid,
        vtkIntArray *CellIds)
{
    vtkIdType NumberOfSeedPoints = SeedGrid->GetNumberOfPoints();
    CellIds->SetNumberOfComponents(1);
    CellIds->SetNumberOfTuples(NumberOfSeedPoints);
    vtkIdType SeedPointCellId;
    double SeedPoint[this->Dimension];
    unsigned long int NumberOfSeedPointsInsideDataGrid = 0;

    // Search Modes
    switch (this->SearchMode)
    {
        // Search Mode Use vtkPointSet
        case SEARCH_MODE_USE_vtkPointSet:
        {
            // Input variables used for vtkPointSet::FindCell
            vtkCell *GuessCell = NULL;   // no guess, global search
            vtkIdType GuessCellId = 0;   // no guess, global search
            double Tolerance = EPSILON;
            // int TempSubId = 0;
            int TempSubId;

            // Output variables of vtkPointSet::FindCell
            double OutputPCoords[3];
            double OutputWeights[4];

            // Locate cells
            for(vtkIdType i=0; i<NumberOfSeedPoints; i++)
            {
                // Get Inquiry point
                SeedGrid->GetPoint(i,SeedPoint);

                // Locate cell
                SeedPointCellId = DataGrid->FindCell(SeedPoint,GuessCell,GuessCellId,Tolerance,TempSubId,OutputPCoords,OutputWeights);

                // Set CellIds in array
                CellIds->SetComponent(i,0,SeedPointCellId);

                // Check points inside the data grid
                if(SeedPointCellId >= 0)
                {
                    NumberOfSeedPointsInsideDataGrid++;
                }

                // Update progress in display
                this->ProgressUpdate(i,NumberOfSeedPoints);
            }

            // end of case I.
            break;
        }

         // Search Mode  Use vtkCellLocator
        case SEARCH_MODE_USE_vtkCellLocator:
        {
            break;
        }

        // Search Mode Use vtkCelltreeLocator
        case SEARCH_MODE_USE_vtkCellTreeLocator:
        {
            break;
        }

        // Search Mode Use vtkHyperOctTree
        case SEARCH_MODE_USE_vtkHyperOctTree:
        {
            break;
        }

        // Search Mode Use vtkModifiedBSPTree
        case SEARCH_MODE_USE_vtkModifiedBSPTree:
        {
            break;
        }

        default:
        {
            vtkErrorMacro("Search mode undefined.");
            return 0;
        }

    } // end of switch

    // Add info to SeedGrid object
    vtkSmartPointer<vtkInformationIntegerKey> NUMBER_OF_SEED_POINTS_INSIDE_DATA_GRID = vtkDataObject::FIELD_NUMBER_OF_TUPLES();
    vtkInformation *SeedGridInfo = SeedGrid->GetInformation();

    if(SeedGridInfo->Has(NUMBER_OF_SEED_POINTS_INSIDE_DATA_GRID))
    {
        SeedGridInfo->Remove(NUMBER_OF_SEED_POINTS_INSIDE_DATA_GRID);
    }

    SeedGridInfo->Set(NUMBER_OF_SEED_POINTS_INSIDE_DATA_GRID,NumberOfSeedPointsInsideDataGrid);

    this->ProgressReset();
    DEBUG(<< this->SearchModeString(this->SearchMode));
    return 1;
}

// ==================
// Initialize Tracers
// ==================

void Seed::InitializeTracers(
        vtkDataSet *DataGrid,
        vtkStructuredPoints *SeedGrid,
        vtkIntArray *CellIds)
{
    // Get Number of Tracers
    vtkSmartPointer<vtkInformation> SeedGridInfo = SeedGrid->GetInformation();
    unsigned long int NumberOfTracers;

    if(SeedGridInfo->Has(vtkDataObject::FIELD_NUMBER_OF_TUPLES()))
    {
        NumberOfTracers = SeedGridInfo->Get(vtkDataObject::FIELD_NUMBER_OF_TUPLES());
    }
    else
    {
        ERROR(<< "NUMBER_OF_SEED_POINTS_INSIDE_DATA_GRID Key is not set.");
        vtkErrorMacro("NUMBER_OF_SEED_POINTS_INSIDE_DATA_GRID Key is not set");
    }

    // Get Number of Seed Grid Points
    unsigned long int NumberOfSeedGridPoints = SeedGrid->GetNumberOfPoints();

    // if(SeedGrid->GetPointData()->HasArray("CellIds") == 0)
    // {
    //     ERROR(<< "CellIds is not in SeedGrid point data.");
    //     vtkErrorMacro(<< "CellIds is not in Seed point data.");
    // }
    // int *pSeedCellIds = vtkIntArray::SafeDownCast(SeedGrid->GetPointData()->GetArray("CellIds"))->GetPointer(0);

    int *pSeedCellIds = CellIds->GetPointer(0);

    // Initialize Tracers
    vtkSmartPointer<vtkDoubleArray> Tracers = vtkSmartPointer<vtkDoubleArray>::New();
    Tracers->SetNumberOfComponents(this->Dimension);
    Tracers->SetNumberOfTuples(NumberOfTracers);
    Tracers->SetName("Tracers");
    // double *pTracers = Tracers->GetPointer(0);
    double *pTracers = Tracers->WritePointer(0,NumberOfTracers * this->Dimension);

    // Initialize TracerIds
    vtkSmartPointer<vtkIntArray> TracerIds = vtkSmartPointer<vtkIntArray>::New();
    TracerIds->SetNumberOfComponents(1);
    TracerIds->SetNumberOfTuples(NumberOfTracers);
    TracerIds->SetName("TracerIds");
    // int *pTracerIds = TracerIds->GetPointer(0);
    int *pTracerIds = TracerIds->WritePointer(0,NumberOfTracers);

    // Get Tracer Points
    unsigned long int TracerCounter = 0;
    for(unsigned long int i=0; i<NumberOfSeedGridPoints; i++)
    {
        if(pSeedCellIds[i] >= 0)
        {
            for(unsigned int j; j<this->Dimension; j++)
            {
                pTracers[this->Dimension * TracerCounter + j] = SeedGrid->GetPoint(i)[j];
            }
            pTracerIds[TracerCounter] = pSeedCellIds[i];
            TracerCounter++;
            this->ProgressUpdate(TracerCounter,NumberOfTracers);
        }
    }

    // Add Number Of Tracers info to DataGrid
    vtkInformation *DataGridInfo = DataGrid->GetInformation();
    if(DataGridInfo->Has(vtkDataObject::FIELD_NUMBER_OF_TUPLES()))
    {
        DataGridInfo->Remove(vtkDataObject::FIELD_NUMBER_OF_TUPLES());
    }
    DataGridInfo->Set(vtkDataObject::FIELD_NUMBER_OF_TUPLES(),NumberOfTracers);

    // Add Tracers to Seed Grid
    SeedGrid->GetFieldData()->AddArray(Tracers);

    // Add TracersId to Seed Grid
    SeedGrid->GetFieldData()->AddArray(TracerIds);

    this->ProgressReset();
    DEBUG(<< "Success");
}

// ==================
// Search Mode String
// ==================

const char * Seed::SearchModeString(SearchModeType SearchMode)
{
    switch(SearchMode)
    {
        // vtkPointSet
        case SEARCH_MODE_USE_vtkPointSet:
        {
            return "SEARCH_MODE_USE_vtkPointSet";
        }

        // vtkCellLocator
        case SEARCH_MODE_USE_vtkCellLocator:
        {
            return "SEARCH_MODE_USE_vtkCellLocator";
        }

        // vtkCellTreeLocator
        case SEARCH_MODE_USE_vtkCellTreeLocator:
        {
            return "SEARCH_MODE_USE_vtkCellTreeLocator";
        }

        // vtkHyperOctTree
        case SEARCH_MODE_USE_vtkHyperOctTree:
        {
            return "SEARCH_MODE_USE_vtkHyperOctTree";
        }

        // vtkModifiedBSPTree
        case SEARCH_MODE_USE_vtkModifiedBSPTree:
        {
            return "SEARCH_MODE_USE_vtkModifiedBSPTree";
        }

        // Not supported type
        default:
        {
            ERROR(<< "Search Mode is not supported.");
            vtkErrorMacro("Search Mode is not supported.");
            return NULL;
        }
    }
}

// ==========================================
// Convert StructuredPoints To StructuredGrid
// ==========================================

void Seed::ConvertStructuredPointsToStructuredGrid(
        vtkStructuredPoints * StructuredPointsData,
        vtkStructuredGrid * StructuredGridData)
{
    // Get Number Of Points
    unsigned int NumberOfPoints = StructuredPointsData->GetNumberOfPoints();

    // Points
    vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
    Points->SetNumberOfPoints(NumberOfPoints);

    // Get Points
    for(unsigned int i=0; i<NumberOfPoints; i++)
    {
        double DoublePoint[3];
        StructuredPointsData->GetPoint(i,DoublePoint);
        Points->SetPoint(i,DoublePoint);
    }

    // Structured Grid
    StructuredGridData->SetPoints(Points);
}
