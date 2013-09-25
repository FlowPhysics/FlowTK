/*
 * =====================================================================================
 *
 *       Filename:  Cache.h
 *
 *    Description:  Cache
 *
 *        Version:  1.0
 *        Created:  09/11/2012 02:14:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Cache_h
#define __Cache_h

// ====================
// Forward Declarations
// ====================

// Complete Declarations
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <BaseFilter.h>

// Imcomplete declarations
class vtkDataSet;
class FilterInformation;

// ===========
// Class Cache
// ===========

class Cache : public vtkMultiBlockDataSetAlgorithm, public BaseFilter
{
    public:
        static Cache * New();
        // vtkTypeRevisionMacro(Cache,vtkMultiBlockDataSetAlgorithm);
        vtkTypeMacro(Cache,vtkMultiBlockDataSetAlgorithm);
        virtual void PrintSelf(ostream &os, vtkIndent indent);
        inline static const char * GetFilterName() { return "Cache"; }

        // Accessors, Mutators
        vtkGetMacro(CacheSize,unsigned int);
        vtkSetMacro(CacheSize,unsigned int);

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
        virtual vtkDataSet * GetInput();
        virtual vtkDataSet * GetInput(int port);
        virtual void SetInput(vtkDataObject *InputDataObject);
        virtual void SetInput(int port, vtkDataObject *InputDataObject);
        virtual void AddInput(vtkDataObject *InputDataObject);
        virtual void AddInput(int port, vtkDataObject *InputDataObject);

    protected:
        Cache();
        virtual ~Cache();

        // Pipeline Executives
        virtual int FillInputPortInformation(int port,vtkInformation *info);
        
        virtual int FillOutputPortInformation(int port,vtkInformation *info);

        virtual int RequestDataObject(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        virtual int RequestInformation(
                vtkInformation *reuqest,
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

        // Member Data
        unsigned int CacheSize;

    private:
        Cache(const Cache &);
        void operator=(const Cache &);
};

#endif
