/*
 * =====================================================================================
 *
 *       Filename:  LCS.h
 *
 *    Description:  Lagrangian Coherent Structures
 *
 *        Version:  1.0
 *        Created:  11/20/2012 06:32:11 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __LCS_h
#define __LCS_h

// ====================
// Forward Declarations
// ====================

// Complete types
#include <vtkDataSetAlgorithm.h>
#include <BaseFilter.h>

// Incomplete types
class FilterInformation;

// =====
// Class
// =====

class LCS : public vtkDataSetAlgorithm, public BaseFilter
{
    public:
        static LCS * New();
        vtkTypeRevisionMacro(LCS,vtkDataSetAlgorithm);
        virtual void PrintSelf(ostream &os, vtkIndent indent);
        virtual inline const char * GetFilterName() { return this->GetStaticFilterName(); }
        static inline const char * GetStaticFilterName() { return "LCS"; }

        // Accessors, Mutators


        // Pipeline Executives
        virtual int ProcessRequest(
                vtkInformation *request,
                vtkInformationVector **input,
                vtkInformationVector *outputVector);

        // Old Pipeline Executives
        virtual vtkStructuredGrid * GetOutput();
        virtual vtkStructuredGrid * GetOutput(int port);
        virtual void SetOutput(vtkDataObject *OutputDataObject);
        virtual void SetOutput(int port, vtkDataObject *OutputDataObject);
        virtual vtkDataSet * GetInput();
        virtual vtkDataSet * GetInput(int port);
        virtual void SetInput(vtkDataObject *InputDataObject);
        virtual void SetInput(int port, vtkDataObject *InputDataObject);
        virtual void AddInput(vtkDataObject *InputDataObject);
        virtual void AddInput(int port, vtkDataObject *InputDataObject);

    protected:
        LCS();
        virtual ~LCS();

        // Pipeeline Executives
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
                vtkInformationVector *outputVector);

        virtual int RequestData(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Member Methods


        // Member Data


    private:
        LCS(const LCS & rhs);
        void operator=(const LCS & rhs);
};

#endif
