/*
 * =====================================================================================
 *
 *       Filename:  FilterInformation.cxx
 *
 *    Description:  Filter Information
 *
 *        Version:  1.0
 *        Created:  05/04/2013 10:58:50 PM
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

#include <FilterInformation.h>

// VTK
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>


// Keys
#include <vtkInformation.h>
#include <vtkInformationKey.h>
#include <vtkInformationDoubleVectorKey.h>
#include <vtkInformationIterator.h>

// ======
// Macros
// ======

vtkStandardNewMacro(FilterInformation);
vtkCxxRevisionMacro(FilterInformation,"$Revision 1.0$");

// ============
//  Constructor
// ============

FilterInformation::FilterInformation()
{
}

// ==========
// Destructor
// ==========

FilterInformation::~FilterInformation()
{
}

// ==========
// Print Self
// ==========

void FilterInformation::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ==============
// TIME_STEPS Key
// ==============

vtkInformationDoubleVectorKey * FilterInformation::TIME_STEPS()
{
    // Key Settings
    const char * TimeStepsKeyName = "TIME_STEPS";
    const char * TimeStepsKeyLocation = FilterInformation::GetStaticFilterName();
    int TimeStepsKeyLength = -1;

    // Create Key
    static vtkInformationDoubleVectorKey * TIME_STEPS_KEY = new vtkInformationDoubleVectorKey(
            TimeStepsKeyName,
            TimeStepsKeyLocation,
            TimeStepsKeyLength);

    return TIME_STEPS_KEY;
}

// ==============
// TIME_RANGE Key
// ==============

vtkInformationDoubleVectorKey * FilterInformation::TIME_RANGE()
{
    // Key Settings
    const char * TimeRangeKeyName = "TIME_RANGE";
    const char * TimeRangeKeyLocation = FilterInformation::GetStaticFilterName();
    int TimeRangeKeyLength = -1;

    // Create Key
    static vtkInformationDoubleVectorKey * TIME_RANGE_KEY = new vtkInformationDoubleVectorKey(
            TimeRangeKeyName,
            TimeRangeKeyLocation,
            TimeRangeKeyLength);

    return TIME_RANGE_KEY;
}

// =====================
// UPDATE_TIME_STEPS Key
// =====================

vtkInformationDoubleVectorKey * FilterInformation::UPDATE_TIME_STEPS()
{
    // Key Settings
    const char * UpdateTimeStepsKeyName = "UPDATE_TIME_STEPS";
    const char * UpdateTimeStepsKeyLocation = FilterInformation::GetStaticFilterName();
    int UpdateTimeStepsKeyLength = -1;

    // Create Key
    static vtkInformationDoubleVectorKey * UPDATE_TIME_STEPS_KEY = new vtkInformationDoubleVectorKey(
            UpdateTimeStepsKeyName,
            UpdateTimeStepsKeyLocation,
            UpdateTimeStepsKeyLength);

    return UPDATE_TIME_STEPS_KEY;
}


// ========
// Get Keys
// ========

// Get Key
vtkInformationKey * FilterInformation::GetKey(
        vtkInformation *InputInformation,
        const char *KeyName)
{
    // Information Iterator
    vtkSmartPointer<vtkInformationIterator> InputInformationIterator = 
        vtkSmartPointer<vtkInformationIterator>::New();
    InputInformationIterator->SetInformationWeak(InputInformation);

    // Iterate over keys
    for(InputInformationIterator->InitTraversal();
            !InputInformationIterator->IsDoneWithTraversal();
            InputInformationIterator->GoToNextItem())
    {
        vtkInformationKey * InputInformationKey = InputInformationIterator->GetCurrentKey();

        if(InputInformationKey == NULL)
        {
            vtkErrorMacro("InputInformationKey is NULL")
        }

        if(!strcmp(InputInformationKey->GetName(),KeyName))
        {
            return InputInformationKey;
        }
    }

    // No key found
    vtkErrorMacro(<< "Key: " << KeyName << " not found.");
    return NULL;
}

// Gey Key
vtkInformationKey * FilterInformation::GetKey(
        vtkInformation *InputInformation,
        vtkInformationKey *SimilarKey)
{
    return this->GetKey(InputInformation,SimilarKey->GetName());
}

// =====================
// Get Double Vector Key
// =====================

// Get Double Vector Key
vtkInformationDoubleVectorKey * FilterInformation::GetDoubleVectorKey(
        vtkInformation *InputInformation,
        const char * KeyName)
{
    // Get Key
    vtkInformationKey * Key = this->GetKey(InputInformation,KeyName);

    // Cast to DoubleKey
    vtkInformationDoubleVectorKey * DoubleVectorKey = 
        vtkInformationDoubleVectorKey::SafeDownCast(Key);

    return DoubleVectorKey;
}

// Get Double Vector Key
vtkInformationDoubleVectorKey * FilterInformation::GetDoubleVectorKey(
        vtkInformation *InputInformation,
        vtkInformationKey *SimilarKey)
{
    return this->GetDoubleVectorKey(InputInformation,SimilarKey->GetName());
}
