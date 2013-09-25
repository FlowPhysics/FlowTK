/*
 * =====================================================================================
 *
 *       Filename:  FlowMap.h
 *
 *    Description:  Flow Map
 *
 *        Version:  1.0
 *        Created:  01/09/2012 02:25:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __FlowMap_h
#define __FlowMap_h

// ====================
// Forward Declarations
// ====================

// Complete Declarations
#include <vtkPolyDataAlgorithm.h>
#include <BaseFilter.h>

// Incomplete Declarations
class vtkStructuredGrid;
class Polynomial;
class RationalNumber;
class FilterInformation;

// =============
// Macros, Enums
// =============

enum IntegratorModeType
{
    INTEGRATOR_MODE_USE_Euler = 0,
    INTEGRATOR_MODE_USE_RungeKutta,
    INTEGRATOR_MODE_USE_AdamsBashforth,
    INTEGRATOR_MODE_NUMBERS
};

// =====
// Class
// =====

class FlowMap : public vtkPolyDataAlgorithm, public BaseFilter
{
    public:
        static FlowMap * New();
        vtkTypeRevisionMacro(FlowMap,vtkPolyDataAlgorithm);
        virtual void PrintSelf(ostream &os, vtkIndent indent);
        inline static const char * GetFilterName() { return "FlowMap"; }

        // Accessors, Mutators
        vtkGetMacro(Dimension,unsigned int);
        vtkSetMacro(Dimension,unsigned int);

        vtkGetMacro(IntegrationTimeStep,double);
        vtkSetMacro(IntegrationTimeStep,double);

        vtkGetMacro(IntegrationDuration,double);
            vtkSetMacro(IntegrationDuration,double);

        vtkGetMacro(IntegratorMode,IntegratorModeType);
        void SetIntegratorMode(IntegratorModeType InputIntegratorMode);

        void SetIntegratorModeToEuler();
        void SetIntegratorModeToRungeKutta();
        void SetIntegratorModeToAdamsBashforth();

        vtkGetMacro(IntegratorOrder,unsigned int);
        vtkSetMacro(IntegratorOrder,unsigned int);

        // Pipeline Executives
        virtual int ProcessRequest(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Old Pipeline Accessors
        virtual vtkStructuredGrid * GetOutput();
        virtual vtkStructuredGrid * GetOutput(int port);
        virtual void SetOutput(vtkDataObject *OutputDataObject);
        virtual void SetOutput(int port, vtkDataObject *OutputDataObject);
        virtual vtkPolyData * GetInput();
        virtual vtkPolyData * GetInput(int port);
        virtual void SetInput(vtkDataObject *InputDataObject);
        virtual void SetInput(int port, vtkDataObject *InputDataObject);
        virtual void AddInput(vtkDataObject *InputDataObject);
        virtual void AddInput(int port, vtkDataObject *InputDataObject);

    protected:
        FlowMap();
        virtual ~FlowMap();

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

        void UpdateInputInformation(
                vtkInformation *inputInfo,
                vtkInformation *outputInfo);

        void UpdateInputDataObject();
        void UpdateInputDataObject(vtkDataObject *InputDataObject);

        // Member Methods
        const char * IntegratorModeString(IntegratorModeType IntegratorMode);
        unsigned int RequiredTimeFrames(IntegratorModeType IntegratorMode);
        void InitializeTimeIndices(vtkInformation *inputInfo0);
        void InitializeIntegratorCoefficients();
        void AdamsBashforthCoefficients(RationalNumber *RationalCoefficients,unsigned int order);
        void RungeKuttaCoefficients(RationalNumber *RationalCoefficients,unsigned int order);
        void InitializeTracers(vtkStructuredGrid *input1);
        double GetCurrentGlobalTime(double CurrentSeedReleaseGlobalTime);
        void AdvectTracers(vtkPolyData *TracerData);

        void PostProcessingTracers(
                vtkPolyData *input0,
                vtkStructuredGrid *input1,
                vtkStructuredGrid *output);

        // Member Data
        unsigned int Dimension;
        IntegratorModeType IntegratorMode;
        unsigned int IntegratorOrder;
        double IntegrationTimeStep;
        double IntegrationDuration;
        double SeedReleaseGlobalTime;
 
        // Internal Member Data
        unsigned long int IntegrationTimeStepIndex;
        unsigned long int IntegrationTimeStepIndexMax;
        double *IntegratorCoefficients;
        vtkPolyData *Tracers;

    private:
        FlowMap(const FlowMap & rhs);
        void operator=(const FlowMap & rhs);
};

#endif
