/*
 * =====================================================================================
 *
 *       Filename:  TestSeed.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/31/2012 11:50:04 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli (sia), 
 *   Organization:  University of California, Berkeley 
 *
 * =====================================================================================
 */

#ifndef __TestSeed_h
#define __TestSeed_h

#include <vtkDataSetAlgorithm.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkStructuredPoints.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

class Test : public vtkDataSetAlgorithm
{
	public:
		static Test * New();
		vtkTypeMacro(Test,vtkDataSetAlgorithm);
		void PrintSelf(ostream os, vtkIndent indent);

	protected:
		Test();
		~Test();

		virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	private:
		Test(const Test &);
		void operator=(const Test &);
};

vtkStandardNewMacro(Test);

void Test::PrintSelf(ostream os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}

Test::Test()
{
}

Test::~Test()
{
}

int Test::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{

	vtkInformation *inputInfo =  inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())) ;

	std::cout << "Inside Test class: " << input->GetNumberOfPoints() << std::endl;
	// outInfo->Print(std::cout);

	return 1;
}

// ====================
// Forward Declarations
// ====================

void CallbackFunction(vtkObject *caller, unsigned long int eventId, void *clientData, void *callData);

#endif
