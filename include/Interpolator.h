/*
 * =====================================================================================
 *
 *       Filename:  Interpolator.h
 *
 *    Description:  Operates computation on data
 *
 *        Version:  1.0
 *        Created:  08/31/2012 04:54:22 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Interpolator_h
#define __Interpolator_h

// =======
// Headers
// =======

// Complete types
#include <vtkDataSetAlgorithm.h>
#include <BaseFilter.h>
#include <vtkInformation.h>

// Incomplete types
class FilterInformation;
class vtkPolyData;
class vtkMultiBlockDataSet;
class vtkIdTypeArray;
class vtkDoubleArray;
class vtkPoints;

// =====
// Types
// =====

enum GridType
{
    UNDEFINED_GRID = 0,
    UNIFORM_GRID,
    HYBRID_GRID,
    NUMBER_OF_GRID_TYPES
};

// ==================
// Class Interpolator
// ==================

class Interpolator : public vtkDataSetAlgorithm, public BaseFilter
{
    public:
        static Interpolator * New();
        vtkTypeRevisionMacro(Interpolator,vtkDataSetAlgorithm);
        virtual void PrintSelf(ostream &os,vtkIndent indent);
        virtual inline const char * GetFilterName() { return this->GetStaticFilterName(); }
        static inline const char * GetStaticFilterName() { return "Interpolator"; }

        // Accessors, Mutators
        vtkGetMacro(SnapToTimeStepTolerance,double);
        vtkSetMacro(SnapToTimeStepTolerance,double);

        // Pipeline Executives
        virtual int ProcessRequest(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Old Pipeline Accessors
        virtual vtkUnstructuredGrid * GetOutput();
        virtual vtkUnstructuredGrid * GetOutput(int port);
        virtual void SetOutput(vtkDataObject *OutputDataObject);
        virtual void SetOutput(int port, vtkDataObject *OutptuDataObject);
        virtual vtkDataSet * GetInput();
        virtual vtkDataSet * GetInput(int port);
        virtual void SetInput(vtkDataObject *InputDataObject);
        virtual void SetInput(int port, vtkDataObject *InputDataObject);
        virtual void AddInput(vtkDataObject *InputDataObject);
        virtual void AddInput(int port, vtkDataObject *InputDataObject);

        // Test
        void Test();

    protected:
        Interpolator();
        virtual ~Interpolator();

        // Pipeline Executives
        virtual int FillInputPortInformation(int port, vtkInformation *info);

        virtual int FillOutputPortInformation(int port, vtkInformation *info);

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
                vtkInformationVector *ouputVector);

        virtual int RequestData(
                vtkInformation *reuqest,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Member Methods
        bool Interpolate(
                vtkMultiBlockDataSet *InputMultiBlockData,
                vtkMultiBlockDataSet *OutputMultiBlockData,
                double *TimeSteps,
                unsigned int TimeStepsLength,
                double *UpdateTimeSteps,
                unsigned int UpdateTimeStepsLength);

        void MultiBlockDataSetToDataObjectArray(
                vtkMultiBlockDataSet *MultiBlockDataSetArray,
                vtkDataObject **DataObjectArray);

        int FindIntervalIndexInArray(
                double InquiryValue,
                double *Array,
                unsigned int ArrayLength);

        double FindTemporalInterpolationCoefficient(
                double UpdateTimeStep,
                double TimeStepLeft,
                double TimeStepRight);

        void LocateTracersInDataGrid(
                vtkDataObject *FirstInputDataObject,
                vtkDataObject **OutputDataObjectsArray,
                unsigned int NumberOfOutputBlocks,
                vtkIdTypeArray **TracersCellIdsArray,
                vtkDoubleArray **TracerParametricCoordinatesInCellArray,
                vtkDoubleArray **TracerWeightsInCellArray);

        void FindGridAdjacencies(
                vtkDataObject *DataObjectGrid,
                vtkIdTypeArray *GridAdjacenciesiArray);

        GridType FindGridTopology(
                vtkDataObject *GridDataObject,
                unsigned int &MaximumNumberOfCellFaces,
                vtkIdTypeArray *NumberOfCellFacesArray);

        vtkIdType FindCellWithLocalSearch(
                double *InquiryPoint,
                unsigned int Dimension,
                vtkIdType GuessCellId,
                vtkDataObject *DataObject,
                double *FoundCellWeights,
                unsigned int NumberOfWeights);

        void SpatioTemporalInterpolation(
                vtkDataObject *InputDataObjectLeft,
                vtkDataObject *InputDataObjectRight,
                vtkDataObject *OutputDataObject,
                vtkIdTypeArray *TracerCellIds,
                vtkDoubleArray *TracerWeightsInCell,
                double TemporalInterpolationCoefficient);

        bool TriangleCellInterpolation(
                const double *InquiryPoint,
                const unsigned int Dimension,
                const double **CellPoints,
                double *InquiryPointLocalCoordinates,
                unsigned int NumberOfLocalCoordinates,
                double *CellPointWeights,
                unsigned int NumberOfWeights);

        bool TetrahedronCellInterpolation(
                const double *InquiryPoint,
                const unsigned int Dimension,
                vtkPoints *CellPoints,
                double *InquiryPointLocalCoordinates,
                unsigned int NumberOfLocalCoordinates,
                double *CellPointWeights,
                unsigned int NumberOfWeights);

        // Member Data
        double SnapToTimeStepTolerance;
        vtkInformation *InterpolatorInfo;

        // Internal Member data
        unsigned int Dimension;
        vtkIdTypeArray *GridAdjacencies;
        GridType DataGridType;

    private:
        Interpolator(const Interpolator &);
        void operator=(const Interpolator &);
};

#endif
