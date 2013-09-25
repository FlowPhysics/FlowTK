/*
 * =====================================================================================
 *
 *       Filename:  Visualization.h
 *
 *    Description:  Visualization (End Of Pipeline)
 *
 *        Version:  1.0
 *        Created:  08/31/2012 05:34:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Visualization_h
#define __Visualization_h

// =======
// Headers
// =======

// Comoplete types
#include <BaseFilter.h>
#include <vtkUnstructuredGridAlgorithm.h>

// Incomplete types
class FilterInformation;

// ==============
// Class Visualization
// ==============

class Visualization : public vtkUnstructuredGridAlgorithm, public BaseFilter
{
	public:
		static Visualization * New();
		vtkTypeRevisionMacro(Visualization,vtkUnstructuredGridAlgorithm);
		virtual void PrintSelf(ostream &os,vtkIndent indent);
		inline static const char * GetFilterName() { return "Visualization"; }

        // Accessors, Mutators

		// Pipeline Executives
		virtual int ProcessRequest(
				vtkInformation *request,
				vtkInformationVector **inputVector,
				vtkInformationVector *outputVector);

	protected:
		Visualization();
		virtual ~Visualization();

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

        virtual int RequestData
            (vtkInformation *request,
             vtkInformationVector **inputVector,
             vtkInformationVector *outputVector);

        // Member Methods
        void VisualizeDataObject(vtkDataObject *DataObject);

        // Member Data
        FilterInformation *VisualizationInformation;

    private:
        Visualization(const Visualization &);
        void operator=(const Visualization &);
};

#endif
