/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  A Complete Pipeline Test
 *
 *        Version:  1.0
 *        Created:  08/31/2012 03:19:07 PM
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

#include <main.h>

// Only for test
#include <vtkExecutive.h>
#include <vtkInformation.h>
#include <vtkStructuredPoints.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkAlgorithm.h>
#include <vtkAlgorithmOutput.h>
#include <vtkDataObject.h>
#include <vtkStructuredGrid.h>

// ====
// Main
// ====

int main(int argc,char *argv[])
{
    // Composite Pipeline //
    vtkSmartPointer<vtkCompositeDataPipeline> prototype = vtkSmartPointer<vtkCompositeDataPipeline>::New();
    vtkAlgorithm::SetDefaultExecutivePrototype(prototype);
    prototype->Delete();

    // Log //
    BaseFilter::SetLogFileName("/home/sia/code/projvtk/Filter/Pipeline/build/log.txt");
    // BaseFilter::LoggerOn();

    // Reader //
    vtkSmartPointer<Reader> reader = vtkSmartPointer<Reader>::New();
    reader->SetNumberOfFileSeries(3);
    reader->SetFileStartNumber(3020);
    reader->SetFileIncrementNumber(20);
    reader->SetFileSeriesPath("/home/sia/data/vtufiles");
    reader->SetFileSeriesBaseName("Patient20Rest-");
    reader->SetFileExtention("vtu");
    reader->DebuggerOn();
    reader->DisplayOn();

    /*
    // Cache
    vtkSmartPointer<Cache> cache = vtkSmartPointer<Cache>::New();
    cache->SetInputConnection(reader->GetOutputPort());
    cache->SetCacheSize(2);
    */

    // Seed //
    vtkSmartPointer<Seed> seed = vtkSmartPointer<Seed>::New();
    seed->SetInputConnection(0,reader->GetOutputPort(0));
    seed->SetDimension(3);
    // seed->SetResolution(40,40,40);
    seed->SetResolution(4,4,4);
    // seed->SetBounds(15,65,105,155,175,225);
    seed->SetBounds(-4,2,4,6,-1,1);
    seed->SetSearchModeToPointSet();
    seed->ProgressOn();
    seed->SetProgressSpeed(1);
    seed->DebuggerOn();
    seed->DisplayOn();
    // seed->Update();

    /*
    // Test Seed I
    vtkSmartPointer<vtkExecutive> SeedExecutive = seed->GetExecutive();
    vtkSmartPointer<vtkInformation> SeedInformation = SeedExecutive->GetOutputInformation(0);

    vtkSmartPointer<vtkStructuredGrid> SeedSP = vtkStructuredGrid::SafeDownCast(SeedInformation->Get(vtkDataObject::DATA_OBJECT()));
    std::cout << BACKGROUND_GREEN << "======= inside main =========" << NONE << std::endl;
    if(SeedSP == NULL)
    {
        std::cout << "SeedSP is NULL" << std::endl;
    }
    std::cout << "Seed in main: number of points: " << SeedSP->GetNumberOfPoints() << std::endl;
    SeedSP->Print(std::cout);
    */

    // /* 
    // Interpolator //
    vtkSmartPointer<Interpolator> interpolator = vtkSmartPointer<Interpolator>::New();
    interpolator->SetInputConnection(0,reader->GetOutputPort(0));
    interpolator->DebuggerOn();
    interpolator->DisplayOn();
    interpolator->ProgressOn();
    interpolator->SetProgressSpeed(1);
    // */

    // /*
    // FlowMap //
    vtkSmartPointer<FlowMap> flowmap = vtkSmartPointer<FlowMap>::New();
    flowmap->SetInputConnection(0,interpolator->GetOutputPort());
    flowmap->SetInputConnection(1,seed->GetOutputPort());
    flowmap->DebuggerOn();
    flowmap->DisplayOn();
    flowmap->ProgressOn();
    flowmap->SetProgressSpeed(1);
    // flowmap->Update();
    // */
    
    /*
    // Deformation //
    vtkSmartPointer<Deformation> deformation = vtkSmartPointer<Deformation>::New();
    deformation->SetInputConnection(flowmap->GetOutputPort());
    deformation->DebuggerOn();
    deformation->DisplayOn();
    // deformation->Update();
    */
    
    /*
    // LCS //
    vtkSmartPointer<LCS> lcs = vtkSmartPointer<LCS>::New();
    lcs->SetInputConnection(deformation->GetOutputPort());
    lcs->DebuggerOn();
    lcs->DisplayOn();
    // lcs->Update();
    */

    // /* 
    // Visualization //
    vtkSmartPointer<Visualization> visualization = vtkSmartPointer<Visualization>::New();
    visualization->SetInputConnection(flowmap->GetOutputPort());
    visualization->DebuggerOn();
    visualization->DisplayOn();
    visualization->Update();
    // */

    /*
    // Derive Pipeline Time
    vtkSmartPointer<vtkStreamingDemandDrivenPipeline> sdd = vtkStreamingDemandDrivenPipeline::SafeDownCast(consumer->GetExecutive());

    unsigned int TimesLength = 2;
    double Times[TimesLength];
    for(unsigned int i=0; i<reader->GetNumberOfFileSeries()-1; i++)
    {
        for(unsigned int j=0; j<TimesLength; j++)
        {
            Times[j] = i+j;
        }
        sdd->SetUpdateTimeSteps(0,Times,2);
        consumer->Update();
        consumer->Print(std::cout);
    }
    */

    /*
    // TEST SEED
    vtkInformation *SeedInfo = seed->GetOutputPortInformation(0);
    if(SeedInfo == NULL)
    {
        std::cout << "SeedInfo is NULL" << std::endl;
    }
    if(! SeedInfo->Has(vtkDataObject::DATA_OBJECT()))
    {
        std::cout << "SeedInfo does not have data object" << std::endl;
    }

    vtkStructuredGrid *SeedOutput = vtkStructuredGrid::SafeDownCast(SeedInfo->Get(vtkDataObject::DATA_OBJECT()));
    std::cout << "End of main, Seed number of points: " << SeedOutput->GetNumberOfPoints() << std::endl;
    */


    HERE_MAIN

    return EXIT_SUCCESS;
}
