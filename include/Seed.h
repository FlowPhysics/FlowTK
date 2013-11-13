/*
 * =====================================================================================
 *
 *       Filename:  Seed.h
 *
 *    Description:  Search for apropriate Seed points
 *
 *        Version:  1.0
 *        Created:  10/16/2012 08:54:27 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Seed_h
#define __Seed_h

// ====================
// Forward Declarations
// ====================

// Complete Declarations
#include <vtkDataSetAlgorithm.h>
#include <BaseFilter.h>
#include <cstdarg>

// Incomplete Declarations
class vtkDataSet;
class vtkIntArray;
class vtkStructuredPoints;

// ============
// Macros/Enums
// ============

enum SearchModeType
{
    SEARCH_MODE_USE_vtkPointSet = 0,
    SEARCH_MODE_USE_vtkCellLocator,
    SEARCH_MODE_USE_vtkCellTreeLocator,
    SEARCH_MODE_USE_vtkHyperOctTree,
    SEARCH_MODE_USE_vtkModifiedBSPTree,
    SEARCH_MODE_NUMBERS
};

// =====
// Class
// =====

class Seed : public vtkDataSetAlgorithm, public BaseFilter
{
    public:
        static Seed * New();
        vtkTypeRevisionMacro(Seed,vtkDataSetAlgorithm);
        virtual void PrintSelf(ostream &os,vtkIndent indent);
        virtual inline const char * GetFilterName() { return this->GetStaticFilterName(); }
        inline static const char * GetStaticFilterName() { return "Seed"; }

        // Accessors, Mutators
        int * GetResolution() const;
        void GetResolution(int *resolution);
        void SetResolution(int *resolution);
        void SetResolution(int resolution0, int resolution1);
        void SetResolution(int resolution0, int reoslution1, int resolution2);

        double * GetBounds() const;
        void GetBounds(double *bounds);
        void SetBounds(double *bounds);
        void SetBounds(double bounds0, double bounds1, double bounds2, double bounds3);
        void SetBounds(double bounds0, double bounds1, double bounds2, double bounds3, double bounds4, double bounds5);

        double * GetSpacing() const;
        void GetSpacing(double *spacing);
        void SetSpacing(double *spacing);
        void SetSpacing(double spacing0, double spacing1);
        void SetSpacing(double spacing0, double spacing1, double spacing2);

        double * GetOrigin() const;
        void GetOrigin(double *origin);
        void SetOrigin(double *origin);
        void SetOrigin(double origin0, double origin1);
        void SetOrigin(double origin0, double origin1, double origin2);

        vtkGetMacro(Dimension,unsigned int);
        vtkSetMacro(Dimension,unsigned int);

        vtkGetMacro(SearchMode,SearchModeType);
        vtkSetMacro(SearchMode,SearchModeType);
        
        void SetSearchModeToPointSet() { this->SearchMode = SEARCH_MODE_USE_vtkPointSet; }
        void SetSearchModeToCellLocator() { this->SearchMode = SEARCH_MODE_USE_vtkCellLocator; }
        void SetSearchModeToCellTreeLocator() { this->SearchMode = SEARCH_MODE_USE_vtkCellTreeLocator; }
        void SetSearchModeToHyperOctTree() { this->SearchMode = SEARCH_MODE_USE_vtkHyperOctTree; }
        void SetSearchModeToModifiedBSPTree() { this->SearchMode = SEARCH_MODE_USE_vtkModifiedBSPTree; }

        // Pipeline Execitives
        virtual int ProcessRequest(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Old Pipeline Accessors
        virtual vtkStructuredPoints * GetOutput();
        virtual vtkStructuredPoints * GetOutput(int port);
        virtual void SetOutput(vtkDataObject *OutputDataObject);
        virtual void SetOutput(int port, vtkDataObject *OutputDataObject);
        virtual vtkDataSet * GetInput();
        virtual vtkDataSet * GetInput(int port);
        virtual void SetInput(vtkDataObject *InputDataObject);
        virtual void SetInput(int port, vtkDataObject *InputDataObject);
        virtual void AddInput(vtkDataObject *InputDataObject);
        virtual void AddInput(int port, vtkDataObject *InputDataObject);

    protected:
        Seed();
        virtual ~Seed();

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
        int SetOutputGridGeometry(
                vtkStructuredPoints *StructuredPoints);

        int FindIntersectionOfGrids(
                vtkDataSet *DataGrid,
                vtkDataSet *SeedGrid,
                vtkIntArray *CellIDs);

        void InitializeTracers(
                vtkDataSet *DataGrid,
                vtkStructuredPoints *SeedGrid,
                vtkIntArray *CellIDs);

        const char * SearchModeString(SearchModeType SearchMode);

        void ConvertStructuredPointsToStructuredGrid(
                vtkStructuredPoints *StructuredPointsData,
                vtkStructuredGrid *StructuredGridData);

        // Member Data
        int *Resolution;
        double *Bounds;
        double *Spacing;
        double *Origin;
        unsigned int Dimension;
        SearchModeType SearchMode;

        // Interal Member data
        double ProcessedTimeStep;

    private:
        Seed(const Seed & rhs);
        void operator=(const Seed & rhs);
};

#endif
