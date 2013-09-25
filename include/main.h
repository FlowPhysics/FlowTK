/*
 * =====================================================================================
 *
 *       Filename:  header.h
 *
 *    Description:  A Complete Pipeline Test
 *
 *        Version:  1.0
 *        Created:  08/31/2012 03:20:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Pipeline_h
#define __Pipeline_h

// =======
// Headers
// =======

// General
#include <vtkSmartPointer.h>
#include <vtkCompositeDataPipeline.h>

// Filters
#include <BaseFilter.h>
#include <Reader.h>
#include <Cache.h>
#include <Seed.h>
#include <Interpolator.h>
#include <FlowMap.h>
#include <Deformation.h>
#include <LCS.h>
#include <Visualization.h>

// ======
// Macros
// ======

// Catch Segmentation Faults
#define HERE_MAIN \
{ \
    std::cout << BACKGROUND_GREEN << "*HERE:  "; \
    std::cout << std::setfill(' ') << std::setw(15) << std::left << "N/A" << " "; \
    std::cout << std::setw(26) << std::left << __FUNCTION__ << " >  "; \
    std::cout << "at line: " << __LINE__; \
    std::cout << NONE << std::endl; \
} \

#endif
