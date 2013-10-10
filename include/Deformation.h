/*
 * =====================================================================================
 *
 *       Filename:  Deformation.h
 *
 *    Description:  Computes Cauchy Green Tensor
 *
 *        Version:  1.0
 *        Created:  11/17/2012 03:36:29 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Deformation_h
#define __Deformation_h

// ====================
// Forward Declarations
// ====================

// Complete
#include <vtkDataSetAlgorithm.h>
#include <BaseFilter.h>

// Incomplete
class FilterInformation;

// =====
// Class
// =====

class Deformation : public vtkDataSetAlgorithm, public BaseFilter
{
    public:
        static Deformation * New();
        // vtkTypeRevisionMacro(Deformation,vtkDataSetAlgorithm);
        vtkTypeMacro(Deformation,vtkDataSetAlgorithm);
        virtual void PrintSelf(ostream &os, vtkIndent indent);
        virtual inline const char * GetFilterName() { return this->GetStaticFilterName(); }
        static inline const char * GetStaticFilterName() { return "Deformation"; }

        // Accessors, Mutators


        // Pipeline Executives
        virtual int ProcessRequest(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Old Pipeeline Accessors
        virtual vtkStructuredGrid * GetOutput();
        virtual vtkStructuredGrid * GetOutput(int port);
        virtual void SetOutput(vtkDataObject *OutputDataObject);
        virtual void SetOutput(int port, vtkDataObject * OutputDataObject);
        virtual vtkDataSet * GetInput();
        virtual vtkDataSet * GetInput(int port);
        virtual void SetInput(vtkDataObject *InputDataObject);
        virtual void SetInput(int port, vtkDataObject *InputDataObject);
        virtual void AddInput(vtkDataObject *InputDataObject);
        virtual void AddInput(int port, vtkDataObject *InputDataObject);

    protected:
        Deformation();
        virtual ~Deformation();

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
                vtkInformationVector *outputVector);

        virtual int RequestData(
            vtkInformation *request,
            vtkInformationVector **inputVector,
            vtkInformationVector *outputVector);

        // Member Methods


        // Member Data


    private:
        Deformation(const Deformation & rhs);
        void operator=(const Deformation & rhs);
};

#endif
