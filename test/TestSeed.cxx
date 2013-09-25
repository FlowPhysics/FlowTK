/*
 * =====================================================================================
 *
 *       Filename:  TestSeed.cxx
 *
 *    Description:  Test Seed class
 *
 *        Version:  1.0
 *        Created:  10/29/2012 12:37:06 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#include <vtkExecutive.h>
#include <vtkInformation.h>
#include <TestSeed.h>

// =======
// Headers
// =======

#include <Seed.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>

// Visualization
#include <vtkVertexGlyphFilter.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindowInteractor.h>

// ====
// Main
// ====

int main(int argc, char *argv[])
{
    // Log
    BaseFilter::SetLogFileName("/home/sia/code/projvtk/Filter/Pipeline/build/log.txt");
//    BaseFilter::LoggerOn();

    // Reader
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
//    reader->SetFileName("/home/sia/data/vtkfiles/restart_vel.2800.vtk");
    reader->SetFileName("/home/sia/data/coronary34_vel.10000.vtk");
    reader->Update();

    // Seed
    vtkSmartPointer<Seed> seed = vtkSmartPointer<Seed>::New();
    seed->SetInputConnection(reader->GetOutputPort());
    seed->SetDimension(3);
    seed->SetResolution(40,40,40);
//    seed->SetBounds(15,65,105,155,175,225);
    seed->SetBounds(-0.5,0.5,0.5,1.5,0.5,2.5);
    seed->SetSearchModeToPointSet();
    seed->ProgressOn();
    seed->SetProgressSpeed(1);
    seed->DebuggerOn();
    seed->DisplayOn();
    seed->Update();

    //////////////////////////////
    // test SEED inside main
    
    vtkExecutive *SeedExecutive = seed->GetExecutive();
    vtkInformation *SeedInfo = SeedExecutive->GetOutputInformation(0);
    vtkStructuredPoints *SeedSP = vtkStructuredPoints::SafeDownCast(SeedInfo->Get(vtkDataObject::DATA_OBJECT()));
    std::cout << "Inside main: " << SeedSP->GetNumberOfPoints() << std::endl;
    // SeedInfo->Print(std::cout);
//    SeedSP->Print(std::cout);

    ///////////////////////////////
    // test Seed in another filter

    vtkSmartPointer<Test> test = vtkSmartPointer<Test>::New();
    test->SetInputConnection(seed->GetOutputPort());
    test->Update();

    //////////////////////////////


    // Write Seed
//    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
//    writer->SetInputConnection(seed->GetOutputPort());
//    writer->SetFileTypeToASCII();
//    writer->SetFileName("/home/sia/Desktop/output.vtk");
//    writer->Update();

    /*
    // Visualization
    vtkSmartPointer<vtkDataSetMapper> DataGridMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    DataGridMapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkSphereSource> GlyphSource = vtkSmartPointer<vtkSphereSource>::New();
    GlyphSource->SetRadius(1);

    vtkSmartPointer<vtkGlyph3D> Glyph = vtkSmartPointer<vtkGlyph3D>::New();
    Glyph->SetSourceConnection(GlyphSource->GetOutputPort());
    Glyph->SetInputConnection(seed->GetOutputPort());
    Glyph->SetColorModeToColorByScalar();
    Glyph->Update();

//    vtkSmartPointer<vtkVertexGlyphFilter> VertexGlyph = vtkSmartPointer<vtkVertexGlyphFilter>::New();
//    VertexGlyph->SetInputConnection(seed->GetOutputPort());

//    vtkSmartPointer<vtkStructuredPointsReader> spreader = vtkSmartPointer<vtkStructuredPointsReader>::New();
//    spreader->SetFileName("/home/sia/Desktop/output.vtk");
//    spreader->Update();
//
//    std::cout << "spreader: " << spreader->GetOutput()->GetNumberOfPoints() << std::endl;

    vtkSmartPointer<vtkDataSetMapper> SeedGridMapper = vtkSmartPointer<vtkDataSetMapper>::New();
//    SeedGridMapper->SetInputConnection(VertexGlyph->GetOutputPort());
    // SeedGridMapper->SetInputConnection(Glyph->GetOutputPort());
    // SeedGridMapper->SetInputConnection(seed->GetOutputPort());
    SeedGridMapper->SetInputConnection(SeedSP->GetProducerPort());

//    vtkStructuredPoints *sp = seed->GetOutput();
    // vtkStructuredPoints *sp = vtkStructuredPoints::SafeDownCast(seed->GetOutputDataObject(0));
//    SeedGridMapper->SetInputConnection(sp->GetProducerPort());
//    SeedGridMapper->SetInputConnection(spreader->GetOutputPort());
//    std::cout << "hahaha: " << sp->GetNumberOfPoints() << std::endl;

    vtkSmartPointer<vtkProperty> DataGridProperty = vtkSmartPointer<vtkProperty>::New();
    DataGridProperty->SetColor(0.5,0.5,1);
    DataGridProperty->SetOpacity(0.7);

    vtkSmartPointer<vtkProperty> SeedGridProperty = vtkSmartPointer<vtkProperty>::New();
    SeedGridProperty->SetColor(1,0,0);
    
    vtkSmartPointer<vtkActor> DataGridActor = vtkSmartPointer<vtkActor>::New();
    DataGridActor->SetMapper(DataGridMapper);
    DataGridActor->SetProperty(DataGridProperty);

    vtkSmartPointer<vtkActor> SeedGridActor = vtkSmartPointer<vtkActor>::New();
    SeedGridActor->SetMapper(SeedGridMapper);
    // SeedGridActor->SetProperty(SeedGridProperty);

    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    // ren->AddActor(DataGridActor);
    ren->AddActor(SeedGridActor);
    // ren->SetBackground(0.9,0.85,1);

    vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();
    renwin->AddRenderer(ren);
    renwin->SetSize(800,600);

    vtkSmartPointer<vtkCallbackCommand> callback = vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction);

    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renwin);
    iren->SetInteractorStyle(style);
    iren->AddObserver(vtkCommand::KeyPressEvent,callback);

    // Start Graphics
    iren->Initialize();
    renwin->Render();
    iren->Start();
    */

    // New Visualization
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(seed->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->AddActor(actor);

    vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();
    renwin->AddRenderer(ren);

    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renwin);

    iren->Initialize();
    renwin->Render();
    iren->Start();


    return EXIT_SUCCESS;
}

// Callback Function
// void CallbackFunction(
//         vtkObject *caller,
//         unsigned long int vtkNotUsed(eventId),
//         void *vtkNotUsed(clientData),
//         void *vtkNotUsed(callData))
// {
//     vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(caller);
//     iren->GetRenderWindow()->Finalize();
//     iren->TerminateApp();
// }
