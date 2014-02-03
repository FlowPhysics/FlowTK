/*
 * =====================================================================================
 *
 *       Filename:  Reader.h
 *
 *    Description:  Reader
 *
 *        Version:  1.0
 *        Created:  08/31/2012 03:21:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Reader_h
#define __Reader_h

// ======
// Macros
// ======

#define ENUM_START 0  // enum start index

// ====================
// Forward Declarations
// ====================

// Complete Types
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <BaseFilter.h>
#include <vtkSmartPointer.h>
#include <map>
#include <vector>

// Incomplete Types
class vtkDataObject;
class File;
class FilterInformation;

// ============
// Class Reader
// ============

class Reader : public vtkMultiBlockDataSetAlgorithm, public BaseFilter
{
    public:
        static Reader * New();
        vtkTypeRevisionMacro(Reader,vtkMultiBlockDataSetAlgorithm);
        virtual void PrintSelf(ostream &os,vtkIndent indent);
        virtual inline const char * GetFilterName() { return "Reader"; }

        // Accessors, Mutators
        vtkGetMacro(NumberOfFileSeries,unsigned int);
        vtkSetMacro(NumberOfFileSeries,unsigned int);

        vtkGetMacro(FileStartNumber,unsigned int);
        vtkSetMacro(FileStartNumber,unsigned int);

        vtkGetMacro(FileIncrementNumber,unsigned int);
        vtkSetMacro(FileIncrementNumber,unsigned int);

        vtkGetStringMacro(FileSeriesPath);
        vtkSetStringMacro(FileSeriesPath);
        
        vtkGetStringMacro(FileSeriesBaseName);
        vtkSetStringMacro(FileSeriesBaseName);

        vtkGetStringMacro(FileExtention);
        vtkSetStringMacro(FileExtention);

        // Pipeline Executives
        virtual int ProcessRequest(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Old Pipeline Accessors
        virtual vtkMultiBlockDataSet * GetOutput();

        virtual vtkMultiBlockDataSet * GetOutput(int port);
        virtual void SetOutput(vtkDataObject *OutputDataObject);
        virtual void SetOutput(int port, vtkDataObject *OutputDataObject);

    protected:
        Reader();
        virtual ~Reader();

        // Pipeline Executives
        virtual int FillInputPortInformation(int port,vtkInformation *info);

        virtual int FillOutputPortInformation(int port,vtkInformation *info);
        
        virtual int RequestDataObject(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        virtual int RequestInformation(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        virtual int RequestUpdateExtent(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        virtual int RequestData(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Member Methods
        int CreateFullFileNames();

        void FindDataTimeSteps();

        bool FindFilesTimeSteps();

        vtkSmartPointer<vtkDataObject> GenericReader(
                const char * FullFileName);

        bool FindMemberIndexInArray(
                double *Array,
                unsigned int ArrayLength,
                double InquiryValue,
                unsigned int &MemberIndex);  // Ouput

        void SortVectorArray(
                double *VectorArray,
                unsigned int VectorArrayLength,
                double *SortedVectorArray);  // Output

        bool FindOutputTimeSteps(
                double *DataTimeSteps,
                unsigned int DataTimeStepsLength,
                double *UpdateTimeSteps,
                unsigned int UpdateTimeStepsLength,
                std::vector<unsigned int> & OutputTimeStepsIndices,  // Output
                std::vector<double> & OutputTimeSteps);  // Output

        // Interface Member Data
        unsigned int NumberOfFileSeries;
        unsigned int FileStartNumber;
        unsigned int FileIncrementNumber;
        char *FileSeriesPath;
        char *FileSeriesBaseName;
        char *FileExtention;

        // Internal Member Data
        FilterInformation *ReaderInformation;
        char **FullFileNames;
        double *DataTimeSteps;
        std::vector<double> OutputTimeStepsVector;
        std::vector<unsigned int> OutputTimeStepsIndicesVector;

    private:
        Reader(const Reader &);
        void operator=(const Reader &);
};

// ==========
// Class File (Container)
// ==========

class File
{
    public:
        // Constructor
        File() {}
        ~File() {}

        // Enum Start index
        // unsigned int ENUM_START;

        // File Type
        enum class FileType : unsigned int
        {
            BINARY = ENUM_START,
            ASCII,
            NUMBER_OF_FILE_TYPES
        };

        // Extension Type
        enum class FileExtensionType
        {
            VTK = ENUM_START,   // Legacy format
            VTI,                // Serial   vtkImageData         (Structured)
            VTP,                // Serial   vtkPolyData          (Structured)
            VTR,                // Serial   vtkRectilinearGrid   (Structured)
            VTS,                // Serial   vtkStructuredGrid    (Structured)
            VTU,                // Serial   vtkUnstructuredGrid  (Unstructured)
            PVTI,               // Parallel vtkImageData         (Structured)
            PVTP,               // Parallel vtkPolyData          (Unstructured)
            PVTR,               // Parallel vtkRectilinearGrid   (Structured)
            PVTS,               // Parallel vtkStructredGrid     (Structured)
            PVTU,               // Parallel vtkUnstructuredgrid  (Unstructured)
            INVALID,            // Not valid format
            NUMBER_OF_FILE_EXTENSIONS
        };

        // DataSet Type
        enum class FileDataSetType
        {
            STRUCTURED_POINTS = ENUM_START,
            STRUCTURED_GRID,
            RECTILINEAR_GRID,
            POLYDATA,
            UNSTRUCTURED_GRID,
            NUMBER_OF_DATASET_TYPES
        };

        // Methods
        static std::map<std::string,File::FileType> InitializeFileTypeMap();
        static std::map<std::string,File::FileExtensionType> InitializeFileExtensionTypeMap();
        static std::map<std::string,File::FileDataSetType> InitializeFileDataSetTypeMap();

        // Member Data (Maps)
        static const std::map<std::string,File::FileType> FileTypeMap;
        static const std::map<std::string,File::FileExtensionType> FileExtensionTypeMap;
        static const std::map<std::string,File::FileDataSetType> FileDataSetTypeMap; 

        // Get File Information
        static File::FileType GetFileType(const char * FullFileName);
        static File::FileExtensionType GetFileExtensionType(const char * FullFileName);
        static File::FileDataSetType GetFileDataSetType(const char * FullFileName);
};

#endif
