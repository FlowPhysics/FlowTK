/*
 * =====================================================================================
 *
 *       Filename:  Reader.cxx
 *
 *    Description:  Reader
 *
 *        Version:  1.0
 *        Created:  08/31/2012 03:38:53 PM
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

// Class
#include <Reader.h>

// For Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// For Keys
#include <FilterInformation.h>
#include <vtkInformationDoubleVectorKey.h>

// General
#include <vtkSmartPointer.h>

// For Reading
#include <vtkXMLFileReadTester.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkGenericDataObjectReader.h>

// Data
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>

// For Strings
#include <sstream>  // stringstream
#include <string>   // string
#include <fstream>  // ifstream
#include <locale>   // isprint
#include <limits>   // numeric_limits
#include <map>      // map

// For DEBUG
#include <vtkInformationRequestKey.h>

// ======
// Macros
// ======

// Time Snap precision
#define EPSILON 1e-4

vtkStandardNewMacro(Reader);
vtkCxxRevisionMacro(Reader,"$Revision 1.0$");

// ===========
// Constructor
// ===========

Reader::Reader()
{
    // Reader Information
    this->ReaderInformation = FilterInformation::New();

    // Member Data
    this->NumberOfFileSeries  = 0;
    this->FileStartNumber     = 0;
    this->FileIncrementNumber = 0;
    this->FileSeriesPath      = NULL;
    this->FileSeriesBaseName  = NULL;
    this->FileExtention       = NULL;
    this->FullFileNames       = NULL;
    this->FileTimeSteps       = NULL;

    // Pipeline
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(2);
}

// ==========
// Destructor
// ==========

Reader::~Reader()
{
    // Reader Information
    if(this->ReaderInformation != NULL)
    {
        this->ReaderInformation->Delete();
    }
    this->ReaderInformation = NULL;

    // FullFileName
    if(this->FullFileNames != NULL)
    {
        delete [] this->FullFileNames;
    }
    this->FullFileNames = NULL;

    // FileTimeSteps
    if(this->FileTimeSteps != NULL)
    {
        delete [] this->FileTimeSteps;
    }
    this->FileTimeSteps = NULL;
}

// ==========
// Print Self
// ==========

void Reader::PrintSelf(ostream &os,vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ==========
// Get Output
// ==========

vtkMultiBlockDataSet * Reader::GetOutput()
{
    return this->GetOutput(0);
}

vtkMultiBlockDataSet * Reader::GetOutput(int port)
{
    return  vtkMultiBlockDataSet::SafeDownCast(this->GetOutputDataObject(port));
}

// ==========
// Set Output
// ==========

void Reader::SetOutput(vtkDataObject *OutputDataObject)
{
    this->SetOutput(0,OutputDataObject);
}

void Reader::SetOutput(int port, vtkDataObject *OutputDataObject)
{
    this->GetExecutive()->SetOutputData(port,OutputDataObject);
}

// ===============
// Process Request
// ===============

int Reader::ProcessRequest(
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

    // Otherwise use supperclass
    return this->Superclass::ProcessRequest(request,inputVector,outputVector);
}

// ===========================
// Fill Input Port Information
// ===========================

int Reader::FillInputPortInformation(int port, vtkInformation *info)
{
    // No input

    DEBUG(<< "Pass");
    return 1;
}

// ============================
// Fill Output Port Information
// ============================

int Reader::FillOutputPortInformation(int port,vtkInformation *info)
{
    // Port 0: vtkMultiBlockDataSet
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkMultiBlockDataSet");
        DEBUG(<< "Success");
        return 1;
    }

    // Port 1: vtkDataObject
    else if(port == 1)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkDataObject");
        DEBUG(<< "Success");
        return 1;
    }

    // Invalid Port
    else
    {
        DEBUG(<< "Failure");
        vtkErrorMacro("Port number is invalid.");
        return 0;
    }
}

// ===================
// Request Data Object
// ===================

int Reader::RequestDataObject(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // Port 0 - vtkMultiBlockDataSet //
 
    // Output
    vtkInformation *outputInfoPort0 = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet *outputPort0 = vtkMultiBlockDataSet::SafeDownCast(outputInfoPort0->Get(vtkDataObject::DATA_OBJECT()));

    // Check info
    if(outputInfoPort0 == NULL)
    {
        ERROR(<< "outputInfoPort0 is NULL.");
        vtkErrorMacro("outputInfoPort1 is NULL.");
        return 0;
    }

    // Create new instance
    if(outputPort0 == NULL)
    {
        outputPort0 = vtkMultiBlockDataSet::New();
        outputInfoPort0->Set(vtkDataObject::DATA_OBJECT(),outputPort0);
        outputPort0->FastDelete();
        this->GetOutputPortInformation(0)->Set(vtkDataObject::DATA_EXTENT_TYPE(),outputPort0->GetExtentType());
    }

    // Port 1 - vtkDataObject //

    // Output
    vtkInformation *outputInfoPort1 = outputVector->GetInformationObject(1);
    vtkDataObject *outputPort1 = vtkDataObject::SafeDownCast(outputInfoPort1->Get(vtkDataObject::DATA_OBJECT()));

    // Check info
    if(outputInfoPort1 == NULL)
    {
        ERROR(<< "outputInfoPort1 is NULL.");
        vtkErrorMacro("outputInfoPort1 is NULL.");
        return 0;
    }

    // Create new instance
    if(outputPort1 == NULL)
    {
        outputPort1 = vtkDataObject::New();
        outputInfoPort1->Set(vtkDataObject::DATA_OBJECT(),outputPort1);
        outputPort1->FastDelete();
        this->GetOutputPortInformation(1)->Set(vtkDataObject::DATA_EXTENT_TYPE(),outputPort1->GetExtentType());
    }

    // Create file series names
    if(!this->FullFileNames)
    {
        this->CreateFullFileNames();
    }

    DEBUG(<< "Success");
    return 1;
}

// ==================
// Request Informaton
// ==================

int Reader::RequestInformation(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // Outputs
    vtkInformation *outputInfoPort0 = outputVector->GetInformationObject(0);
    vtkInformation *outputInfoPort1 = outputVector->GetInformationObject(1);

    // Check outputInfo is not empty
    if(outputInfoPort0 == NULL)
    {
        vtkErrorMacro("outputInfoPort0 is NULL.");
        return 0;
    }

    if(outputInfoPort1 == NULL)
    {
        vtkErrorMacro("outputInfoPort1 is NULL.");
        return 0;
    }

    // 1- DATA TIME STEPS //

    // Check if Key Exists in outputInfo
    bool InitializeDataTimeStepsKey = false;

    // Check if key is included in outputInfo
    if(!outputInfoPort0->Has(FilterInformation::DATA_TIME_STEPS()))
    {
        InitializeDataTimeStepsKey = true;
    }
    else
    {
        // Check if key length is correct
        if(outputInfoPort0->Length(FilterInformation::DATA_TIME_STEPS()) != 
                static_cast<int>(this->NumberOfFileSeries))
        {
            InitializeDataTimeStepsKey = true;
        }

        // Check if key values are not empty
        if(outputInfoPort0->Get(FilterInformation::DATA_TIME_STEPS()) == NULL)
        {
            InitializeDataTimeStepsKey = true;
        }
    }

    // Initialize DataTimeSteps Key
    if(InitializeDataTimeStepsKey == true)
    {
        // DataTimeSteps Length
        int DataTimeStepsLength = this->NumberOfFileSeries;

        // DataTimeSteps Values
        double DataTimeSteps[DataTimeStepsLength];

        // Check if files have default data time steps
        bool DataTimeStepsAvailable = this->GetFilesTimeSteps();

        if(DataTimeStepsAvailable == true)
        {
            // Use default time steps in files
            for(unsigned int i=0; static_cast<int>(i)<DataTimeStepsLength; i++)
            {
                DataTimeSteps[i] = this->FileTimeSteps[i];
            }
        }
        else
        {
            // Use file index for time steps
            for(unsigned int i=0; static_cast<int>(i)<DataTimeStepsLength; i++)
            {
                DataTimeSteps[i] = double(i);
            }
        }

        // Add values to key in outputInfo
        outputInfoPort0->Set(FilterInformation::DATA_TIME_STEPS(),DataTimeSteps,DataTimeStepsLength);
        outputInfoPort1->Set(FilterInformation::DATA_TIME_STEPS(),DataTimeSteps,DataTimeStepsLength);

        // Debug
        DISPLAY(DataTimeSteps,this->NumberOfFileSeries);
    }

    // 2- DATA TIME RANGE //

    // Check if Key exists
    bool InitializeDataTimeRangeKey = false;

    // Check if key is included in outputInfo
    if(!outputInfoPort0->Has(FilterInformation::DATA_TIME_RANGE()))
    {
        InitializeDataTimeRangeKey = true;
    }
    else
    {
        // Check if key length is correct
        if(outputInfoPort0->Length(FilterInformation::DATA_TIME_RANGE()) != 2)
        {
            InitializeDataTimeRangeKey = true;
        }

        // Check if key values are empty
        if(outputInfoPort0->Get(FilterInformation::DATA_TIME_RANGE()) == NULL)
        {
            InitializeDataTimeRangeKey = true;
        }
    }

    // Initialize DataTimeRange Key
    if(InitializeDataTimeRangeKey == true)
    {
        // Get TimeSteps
        double *DataTimeSteps = outputInfoPort0->Get(FilterInformation::DATA_TIME_STEPS());
        int DataTimeStepsLength = outputInfoPort0->Length(FilterInformation::DATA_TIME_STEPS());

        // Create TimeRange values
        double DataTimeRange[2];
        DataTimeRange[0] = DataTimeSteps[0];
        DataTimeRange[1] = DataTimeSteps[DataTimeStepsLength-1];

        // Add values to key in outputInfo
        outputInfoPort0->Set(FilterInformation::DATA_TIME_RANGE(),DataTimeRange,2);
        outputInfoPort1->Set(FilterInformation::DATA_TIME_RANGE(),DataTimeRange,2);

        // Debug
        DISPLAY(DataTimeRange,2);
    }

    // Debug //
    DEBUG(<< "Success");

    return 1;
}

// =====================
// Request Update Extent
// =====================

int Reader::RequestUpdateExtent(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *vtkNotUsed(outputVector))
{
    // No input

    DEBUG(<< "Pass");
    return 1;
}

// ============
// Request Data
// ============

int Reader::RequestData(
        vtkInformation *request,
        vtkInformationVector **vtkNotUsed(inputVector),
        vtkInformationVector *outputVector)
{
    // output Port 0

    vtkInformation *outputInfoPort0 = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet *outputPort0 = 
        vtkMultiBlockDataSet::SafeDownCast(outputInfoPort0->Get(vtkDataObject::DATA_OBJECT()));

    // Check info is not NULL
    if(outputInfoPort0 == NULL)
    {
        ERROR(<< "outputInfo of port 0 is NULL.");
        vtkErrorMacro("outputInfo of port 0 is NULL.");
        return 0;
    }

    // Check output port 0 is not NULL
    if(outputPort0 == NULL)
    {
        ERROR(<< "output of port 0 is NULL.");
        vtkErrorMacro("output of port 0 is NULL.");
        return 0;
    }

    // Check outputInfo of port 0 has TIME_STEPS key
    if(!outputInfoPort0->Has(FilterInformation::TIME_STEPS()))
    {
        ERROR(<< "outputInfo of port 0 does not have TIME_STEPS key.");
        vtkErrorMacro("outputInfo of port 0 does not have TIME_STEPS key.");
        return 0;
    }

    // Data Time Steps
    double *DataTimeSteps = outputInfoPort0->Get(FilterInformation::DATA_TIME_STEPS());
    unsigned int DataTimeStepsLength = outputInfoPort0->Length(FilterInformation::DATA_TIME_STEPS());

    // Check Data Time Steps
    if(DataTimeSteps == NULL)
    {
        ERROR(<< "DataTimeSteps is NULL.");
        vtkErrorMacro("DataTimeSteps is NULL.");
        return 0;
    }

    if(DataTimeStepsLength < 1)
    {
        ERROR(<< "DataTimeStepsLength is zero.");
        vtkErrorMacro("DataTimeStepsLength is zero.");
        return 0;
    }

    // Check outputInfo of port 0 has UPDATE_TIME_STEPS key
    if(!outputInfoPort0->Has(FilterInformation::UPDATE_TIME_STEPS()))
    {
        ERROR(<< "outputInfo of port 0 does not have UPDATE_TIME_STEPS key.");
        vtkErrorMacro("outputInfoPort0 does not have UPDATE_TIME_STEPS key.");
        return 0;
    }

    // Update Time Steps
    double *UpdateTimeSteps = outputInfoPort0->Get(FilterInformation::UPDATE_TIME_STEPS());
    unsigned int UpdateTimeStepsLength = outputInfoPort0->Length(FilterInformation::UPDATE_TIME_STEPS());

    // Check Update Time Steps
    if(UpdateTimeSteps == NULL)
    {
        ERROR(<< "UpdateTimeSteps is NULL.");
        vtkErrorMacro("UpdateTimeSteps is NULL.");
        return 0;
    }

    if(UpdateTimeStepsLength < 1)
    {
        ERROR(<< "UpdateTimeStepsLength is zero");
        vtkErrorMacro("Update Time Steps Length should be at least 1.");
        return 0;
    }

    // Output Block DataSet Length
    outputPort0->SetNumberOfBlocks(UpdateTimeStepsLength);

    // Read Data //

    // Iterate over Update time steps
    for(unsigned int UpdateTimeStepsIterator = 0;
            UpdateTimeStepsIterator < UpdateTimeStepsLength;
            UpdateTimeStepsIterator++)
    {
        // Index of Updated Time Step
        unsigned int FileIndex;
        bool UniqueIndexFound = this->FindIndexInVectorArray(
                DataTimeSteps,
                DataTimeStepsLength,
                UpdateTimeSteps[UpdateTimeStepsIterator],
                FileIndex);                // Output

        // Check if unique index is found
        if(UniqueIndexFound == false)
        {
            ERROR(<< "Unique index not found.");
            vtkErrorMacro("Reader >> Unique index not found!");
            return 0;
        }

        // Generic Reader
        vtkSmartPointer<vtkDataObject> outputDataObject = this->GenericReader(this->FullFileNames[FileIndex]);
        // vtkDataObject * outputDataObject = this->GenericReader(this->FullFileNames[FileIndex]); 
        // Use the above line instead of two lines above, WITH uncommenting REGISTER(0) in GenericReader method.
        // and removing the vtkSmartPointer in SIGNATURE of GenericReader method.
        
        // Check valid data
        if(outputDataObject == NULL)
        {
            ERROR(<< "outputDataObject is NULL.");
            vtkErrorMacro("outputDataObject is NULL.");
            return 0;
        }
        else
        {
            // Check non-zero number of points
            vtkSmartPointer<vtkDataSet> outputDataSet = vtkDataSet::SafeDownCast(outputDataObject);
            if(outputDataSet->GetNumberOfPoints() < 1)
            {
                ERROR(<< "outputDataObject has no points.");
                vtkErrorMacro("outputDataObject has no points.");
                return 0;
            }
        }

        DEBUG(<< "Read: " << this->FullFileNames[FileIndex])

        // Add a leaf
        outputPort0->SetBlock(UpdateTimeStepsIterator,outputDataObject);
    }

    // Add time information
    outputPort0->GetInformation()->Set(FilterInformation::UPDATE_TIME_STEPS(),UpdateTimeSteps,UpdateTimeStepsLength);

    /// Output Port 1 ///
    
    // vtkInformation *outputInfoPort1 = outputVector->GetInformationObject(1);
    // vtkDataSet *outputPort1 = vtkDataSet::SafeDownCast(outputInfoPort1->Get(vtkDataObject::DATA_OBJECT()));

    // outputPort1 =  // TODO

    // Set Output TIME STEPS //

    // Put UpdateTimeSteps to OutputUpdateTimeSteps
    double *OutputTimeSteps = UpdateTimeSteps;
    unsigned int OutputTimeStepsLength = UpdateTimeStepsLength;
    outputInfoPort0->Set(FilterInformation::TIME_STEPS(),OutputTimeSteps,OutputTimeStepsLength);

    // Debug
    DISPLAY(OutputTimeSteps,OutputTimeStepsLength);


    // Debug //
    DISPLAY(DataTimeSteps,DataTimeStepsLength);
    DISPLAY(UpdateTimeSteps,UpdateTimeStepsLength);
    DEBUG(<< "Success");

    return 1;
}

// ======================
// Create Full File Names
// ======================

int Reader::CreateFullFileNames()
{
    if(this->NumberOfFileSeries == 0) 
    {
        vtkErrorMacro("Reader >> NumberOfFileSeries should be greater than zero.");
        return 0;
    }
    
    this->FullFileNames = new char *[this->NumberOfFileSeries];

    unsigned int StrLength = strlen(this->FileSeriesPath) + strlen(this->FileSeriesBaseName) + strlen(this->FileExtention) + 
        log10(this->FileStartNumber + this->FileIncrementNumber * (this->NumberOfFileSeries-1)) + 3;

    for(unsigned int i=0, FileNumber = this->FileStartNumber; i<this->NumberOfFileSeries; i++, FileNumber +=this->FileIncrementNumber)
    {
        std::stringstream StrStream;
        StrStream << this->FileSeriesPath << "/" << this->FileSeriesBaseName << FileNumber << "." << this->FileExtention;
        FullFileNames[i] = new char[StrLength+1];
        strcpy(FullFileNames[i],StrStream.str().c_str());
        StrStream.str("");
    }

    DEBUG(<< "Success");
    return 1;
}

// ====================
// Get Files Time Steps
// ====================

// Description
// Returns false if time sterps are not found or not in valid format

bool Reader::GetFilesTimeSteps()
{
    // Create Full File Names
    if(this->FullFileNames == NULL)
    {
        this->CreateFullFileNames();
    }

    // Iterate over all files
    for(unsigned int i=0; i<this->NumberOfFileSeries; i++)
    {
        // Open File
        std::ifstream File;
        File.open(this->FullFileNames[0],std::ios::in);
        if(!File.is_open())
        {
            vtkErrorMacro("Can not open the file.");
        }

        // Read File (Line by line)
        File.seekg(0,std::ios::beg);
        std::string Line = "";
        File.ignore(std::numeric_limits<std::streamsize>::max(),'\n');  // Skip first line till \n
        std::getline(File,Line);  // Read second line


        // Check if Time exists
        bool ValidTimeFormat = true;
        std::string ValidNumerics = "0123456789";
        size_t FoundLastNumber = Line.find_last_of(ValidNumerics);
        if(FoundLastNumber == std::string::npos)
        {
            // No number found
            return false;
        }
        else
        {
            // Alphabets after last number
            std::string AlphabetsAfterLastNumber = Line.substr(FoundLastNumber+1,std::string::npos);
            if(AlphabetsAfterLastNumber.length() > 1)
            {
                ValidTimeFormat = false;
            }
            else if(AlphabetsAfterLastNumber.length() == 1)
            {
                // Find if there is alphabet after numbers, excluding dfs.,;
                // s for sec, d for integrs, f for floats are allowed at the end
                size_t FoundLastAlphabet = AlphabetsAfterLastNumber.find_last_not_of(".,;sfd");
                if(FoundLastAlphabet != std::string::npos)
                {
                    ValidTimeFormat = false;
                }
            }
        }

        // Check if time is a valid format
        if(ValidTimeFormat == false)
        {
            return false;
        }

        // Extract last numbers
        std::string ValidNumbers = ValidNumerics + ".-+";
        size_t FoundAlphabetBeforeLastNumber = Line.find_last_not_of(ValidNumbers);
        std::string FileTimeAsString = Line.substr(FoundAlphabetBeforeLastNumber+1,FoundLastNumber);

        // String to double conversion
        double FileTime = atof(FileTimeAsString.c_str());

        // Store in Time Array
        if(this->FileTimeSteps == NULL)
        {
            this->FileTimeSteps = new double[this->NumberOfFileSeries];
        }

        this->FileTimeSteps[i] = FileTime;

        // Close File
        File.close();

        HERE
    }

    HERE

    return true;
}

// ==========================
// Find Index In Vector Array
// ==========================

// Description:
// Find index of a member of an array.
// If ArrayMember is not in range of array or not close to array members it
// returns false, otherwise it returns true.

bool Reader::FindIndexInVectorArray(
        double * VectorArray,
        unsigned int VectorArrayLength,
        double ArrayMember,
        unsigned int &MemberIndex)
{
    // Count How many index will be found
    unsigned int IndexFound = 0;

    // Iterate over VectorArray
    for(unsigned int i=0; i<VectorArrayLength; i++)
    {
        if(abs(VectorArray[i]-ArrayMember) < EPSILON)
        {
            IndexFound++;
            MemberIndex = i;
        }
    }

    // Check unique index
    if(IndexFound != 1)
    {
        return false;
    }
    else
    {
        return true;
    }
}

// ==============
// Generic Reader
// ==============

vtkSmartPointer<vtkDataObject> Reader::GenericReader(const char *FullFileName)
{
    // Determine VTK or XML file
    vtkXMLFileReadTester *XMLFileReadTester = vtkXMLFileReadTester::New();
    XMLFileReadTester->SetFileName(FullFileName);

    int FileIsXML = XMLFileReadTester->TestReadFile(); 
    
    // Freed memory
    XMLFileReadTester->Delete();

    if(FileIsXML == 1)
    {
        // File is XML
        vtkSmartPointer<vtkXMLGenericDataObjectReader> XMLReader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
        XMLReader->SetFileName(FullFileName);
        XMLReader->Update();

        vtkSmartPointer<vtkDataObject> DataObject = XMLReader->GetOutput();

        if(DataObject == NULL)
        {
            vtkErrorMacro("DataObject is NULL.");
            return NULL;
        }

        // DataObject->Register(0);

        return DataObject;
    }
    else
    {
        // File is not XML
        vtkSmartPointer<vtkGenericDataObjectReader> Reader = 
            vtkSmartPointer<vtkGenericDataObjectReader>::New();
        Reader->SetFileName(FullFileName);
        Reader->Update();

        vtkSmartPointer<vtkDataObject> DataObject = Reader->GetOutput();

        if(DataObject == NULL)
        {
            vtkErrorMacro("DataObject is NULL.");
            return NULL;
        }

        // DataObject->Register(0);

        return DataObject;
    }
}

// =============
// Get File Type
// =============

// Determines if the first input file is BINARy or ASCII.

File::FileType File::GetFileType(const char * FullFileName)
{
    // Open file
    std::ifstream FirstFile;
    FirstFile.open(FullFileName,std::ios::binary);
    if(!FirstFile.is_open())
    {
        std::cerr << "Can not read the file." << std::endl;
    }

    // Read file (byte by byte)
    FirstFile.seekg(0,std::ios::beg);
    size_t NumberOfBytesToRead = 2000;
    char Array[NumberOfBytesToRead];
    FirstFile.read(Array,NumberOfBytesToRead);  // Read by bytes

    // Look for Printable characters
    for(unsigned int i=0; i<NumberOfBytesToRead; i++)
    {
        // Exclude "new line" character (it is not printable)
        if(!std::isprint(Array[i]) && (Array[i] != '\n'))
        {
            return File::FileType::BINARY;
        }
    }

    FirstFile.close();
    return File::FileType::ASCII;
}

// =======================
// Get File Extention Type
// =======================

File::FileExtensionType File::GetFileExtensionType(const char *FullFileName)
{
    // Filename as String
    std::string FullFileNameAsString("FullFileName");

    // Find last dot
    size_t FoundLastDot = FullFileNameAsString.find_last_of(".");

    // Check valid dot position
    if(FoundLastDot == std::string::npos)
    {
        std::cerr << "Filename does not have extenssion." << std::endl;
        return File::FileExtensionType::INVALID;
    }

    // Extract string after dot
    std::string FileExtensionAsString = FullFileNameAsString.substr(FoundLastDot+1,std::string::npos);

    // Determine Extension Type
    std::map<std::string,File::FileExtensionType>::const_iterator FileExtensionTypeMapIterator;
    FileExtensionTypeMapIterator = File::FileExtensionTypeMap.find(FileExtensionAsString);

    // Check if Extension is valid
    if(FileExtensionTypeMapIterator == File::FileExtensionTypeMap.end())
    {
        std::cerr << "File extenstion type: " << FileExtensionAsString << " is invalid." << std::endl;
        return File::FileExtensionType::INVALID;
    }
    else
    {
        return FileExtensionTypeMapIterator->second;
    }
}

// =====================
// Get File DataSet Type
// =====================

File::FileDataSetType File::GetFileDataSetType(const char * FullFileName)
{
    // TODO
    return File::FileDataSetType::UNSTRUCTURED_GRID;
}

// ======================
// Initialize static maps
// ======================

const std::map<std::string,File::FileType> File::FileTypeMap
      = File::InitializeFileTypeMap();

const std::map<std::string,File::FileExtensionType> File::FileExtensionTypeMap
      = File::InitializeFileExtensionTypeMap();

const std::map<std::string,File::FileDataSetType> File::FileDataSetTypeMap
      = File::InitializeFileDataSetTypeMap();

// ========================
// Initialize File Type Map
// ========================

std::map<std::string,File::FileType> File::InitializeFileTypeMap()
{
    // Create map
    std::map<std::string,File::FileType> FileTypeMap;

    // Insert Pairs
    FileTypeMap["BINARY"] = File::FileType::BINARY;
    FileTypeMap["ASCII" ] = File::FileType::ASCII;
    
    return FileTypeMap;
}

// ==================================
// Initialize File Extension Type Map
// ==================================

std::map<std::string,File::FileExtensionType> File::InitializeFileExtensionTypeMap()
{
    // Create map
    std::map<std::string,File::FileExtensionType> FileExtensionTypeMap;
    
    // Insert Pairs
    FileExtensionTypeMap["VTK" ] = File::FileExtensionType::VTK;
    FileExtensionTypeMap["VTI" ] = File::FileExtensionType::VTI;
    FileExtensionTypeMap["VTP" ] = File::FileExtensionType::VTP;
    FileExtensionTypeMap["VTR" ] = File::FileExtensionType::VTR;
    FileExtensionTypeMap["VTS" ] = File::FileExtensionType::VTS;
    FileExtensionTypeMap["VTU" ] = File::FileExtensionType::VTU;
    FileExtensionTypeMap["PVTI"] = File::FileExtensionType::PVTI;
    FileExtensionTypeMap["PVTP"] = File::FileExtensionType::PVTP;
    FileExtensionTypeMap["PVTR"] = File::FileExtensionType::PVTR;
    FileExtensionTypeMap["PVTS"] = File::FileExtensionType::PVTS;
    FileExtensionTypeMap["PVTU"] = File::FileExtensionType::PVTU;

    return FileExtensionTypeMap;
}

// ================================
// Initialize File DataSet Type Map
// ================================

std::map<std::string,File::FileDataSetType> File::InitializeFileDataSetTypeMap()
{
    // Create Map
    std::map<std::string,File::FileDataSetType> FileDataSetTypeMap;

    // Insert Pairs
    FileDataSetTypeMap["STRUCTURED_POINTS"] = File::FileDataSetType::STRUCTURED_POINTS;
    FileDataSetTypeMap["STRUCTURED_GRID"  ] = File::FileDataSetType::STRUCTURED_GRID;
    FileDataSetTypeMap["RECTILINEAR_GRID" ] = File::FileDataSetType::RECTILINEAR_GRID;
    FileDataSetTypeMap["POLYDATA"         ] = File::FileDataSetType::POLYDATA;
    FileDataSetTypeMap["UNSTRUCTURED_GRID"] = File::FileDataSetType::UNSTRUCTURED_GRID;

    return FileDataSetTypeMap;
}
