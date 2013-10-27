/*
 * =====================================================================================
 *
 *       Filename:  FilterInformation.h
 *
 *    Description:  Filter Information
 *
 *        Version:  1.0
 *        Created:  05/04/2013 10:59:18 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __FilterInformation_h
#define __FilterInformation_h

// =======
// Headers
// =======

// Complete Types
#include <vtkInformation.h>
#include <BaseFilter.h>

// Incomplete Types
class vtkInformationKey;
class vtkInformationDoubleVectorKey;

// =================
// FilterInformation
// =================

class FilterInformation : public vtkInformation, public BaseFilter
{
    public:
        static FilterInformation * New();
        vtkTypeRevisionMacro(FilterInformation,vtkInformation);
        virtual void PrintSelf(ostream &os, vtkIndent Indent);
        virtual inline const char * GetFilterName() { return this->GetStaticFilterName(); }
        static const char * GetStaticFilterName() { return "FilterInformation"; }

        // Member Methods
        vtkInformationKey * GetKey(
                vtkInformation *InputInformation,
                const char *KeyName);

        vtkInformationKey * GetKey(
                vtkInformation *InputInformation,
                vtkInformationKey *SimilarKey);

        vtkInformationDoubleVectorKey * GetDoubleVectorKey(
                vtkInformation * InputInformation,
                const char * KeyName);

        vtkInformationDoubleVectorKey * GetDoubleVectorKey(
                vtkInformation *InputInformation,
                vtkInformationKey *SimilarKey);

		// Keys
		static vtkInformationDoubleVectorKey * TIME_STEPS();
		static vtkInformationDoubleVectorKey * TIME_RANGE();
		static vtkInformationDoubleVectorKey * DATA_TIME_STEPS();
		static vtkInformationDoubleVectorKey * UPDATE_TIME_STEPS();
        
    protected:
        FilterInformation();
        virtual ~FilterInformation();

        // Member Data
        vtkInformation *Information;

    private:
        FilterInformation(const FilterInformation &rhs);
        void operator=(const FilterInformation &rhs);
};

#endif

