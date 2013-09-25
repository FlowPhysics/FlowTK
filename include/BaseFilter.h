/*
 * =====================================================================================
 *
 *       Filename:  BaseFilter.hpp
 *
 *    Description:  Base Filter
 *
 *        Version:  1.0
 *        Created:  09/12/2012 02:55:38 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __BaseFilter_hpp
#define __BaseFilter_hpp

// =======
// Headers
// =======

#include <vtkSetGet.h>   // vtk Macros
#include <iomanip>       // Formatting output
#include <fstream>       // write log file
#include <string>

// ======
// Macros
// ======

// Detect OS
#if defined(__WIN32__) || defined(__WIN32) || defined(_WIN32) || defined(WIN32)
    #define __WINDOWS__
#endif


// Colors
#if defined(__WINDOWS__)                   // No ANSI code support in windows
    const char FOREGROUND_BLACK[]          = "";
    const char FOREGROUND_RED[]            = "";
    const char FOREGROUND_GREEN[]          = "";
    const char FOREGROUND_YELLOW[]         = "";
    const char FOREGROUND_BLUE[]           = "";
    const char FOREGROUND_MAGENTA[]        = "";
    const char FOREGROUND_CYAN[]           = "";
    const char FOREGROUND_WHITE[]          = "";

    const char FOREGROUND_LIGHT_BLACK[]    = "";
    const char FOREGROUND_LIGHT_RED[]      = "";
    const char FOREGROUND_LIGHT_GREEN[]    = "";
    const char FOREGROUND_LIGHT_YELLOW[]   = "";
    const char FOREGROUND_LIGHT_BLUE[]     = "";
    const char FOREGROUND_LIGHT_MAGENTA[]  = "";
    const char FOREGROUND_LIGHT_CYAN[]     = "";
    const char FOREGROUND_LIGHT_WHITE[]    = "";

    const char BACKGROUND_BLACK[]          = "";
    const char BACKGROUNG_RED[]            = "";
    const char BACKGROUND_GREEN[]          = "";
    const char BACKGROUND_YELLOW[]         = "";
    const char BACKGROUND_BLUE[]           = "";
    const char BACKGROUND_MAGENTA[]        = "";
    const char BACKGROUND_CYAN[]           = "";
    const char BACKGROUND_WHITE[]          = "";

    const char NONE[]                      = "";
#else
    const char FOREGROUND_BLACK[]          = "\033[0;30m";
    const char FOREGROUND_RED[]            = "\033[0;31m";
    const char FOREGROUND_GREEN[]          = "\033[0;32m";
    const char FOREGROUND_YELLOW[]         = "\033[0;33m";
    const char FOREGROUND_BLUE[]           = "\033[0;34m";
    const char FOREGROUND_MAGENTA[]        = "\033[0;35m";
    const char FOREGROUND_CYAN[]           = "\033[0;36m";
    const char FOREGROUND_WHITE[]          = "\033[0;37m";

    const char FOREGROUND_LIGHT_BLACK[]    = "\033[1;30m";
    const char FOREGROUND_LIGHT_RED[]      = "\033[1;31m";
    const char FOREGROUND_LIGHT_GREEN[]    = "\033[1;32m";
    const char FOREGROUND_LIGHT_YELLOW[]   = "\033[1;33m";
    const char FOREGROUND_LIGHT_BLUE[]     = "\033[1;34m";
    const char FOREGROUND_LIGHT_MAGENTA[]  = "\033[1;35m";
    const char FOREGROUND_LIGHT_CYAN[]     = "\033[1;36m";
    const char FOREGROUND_LIGHT_WHITE[]    = "\033[1;37m";

    const char BACKGROUND_BLACK[]          = "\033[0;40m";
    const char BACKGROUND_RED[]            = "\033[0;41m";
    const char BACKGROUND_GREEN[]          = "\033[0;42m";
    const char BACKGROUND_YELLOW[]         = "\033[0;43m";
    const char BACKGROUND_BLUE[]           = "\033[0;44m";
    const char BACKGROUND_MAGENTA[]        = "\033[0;45m";
    const char BACKGROUND_CYAN[]           = "\033[0;46m";
    const char BACKGROUND_WHITE[]          = "\033[0;47m";

    const char NONE[]                      = "\033[0m"   ;
#endif


// Output Stream
#define OutputStream *BaseFilter::BaseFilterOutputStream


// DEBUG
#define DEBUG(...) \
{ \
    if(this->DebugStatus) \
    { \
        OutputStream << "DEBUG:  "; \
        OutputStream << std::setfill(' ') << std::setw(15) << std::left << this->GetClassName() << " "; \
        OutputStream << std::setw(26) << std::left << __FUNCTION__ << " >  "; \
        if(strstr(#__VA_ARGS__,"\"Success\"") != NULL) \
        { \
            OutputStream << FOREGROUND_GREEN; \
        } \
        else if(strstr(#__VA_ARGS__,"\"Failure\"") != NULL) \
        { \
            OutputStream << FOREGROUND_RED; \
        } \
        else if(strstr(#__VA_ARGS__,"\"Pass\"") != NULL) \
        { \
            OutputStream << FOREGROUND_YELLOW; \
        } \
        OutputStream __VA_ARGS__; \
        OutputStream << NONE << std::endl; \
    } \
} \


// DISPLAY
#define DISPLAY(FirstArg,...) \
{ \
    if(this->DisplayStatus) \
    { \
        DisplayFunction(this->GetClassName(),__FUNCTION__,#FirstArg,(FirstArg),##__VA_ARGS__); \
    } \
} \


// PRINT
#define PRINT(...) \
{ \
    OutputStream << FOREGROUND_BLUE << "PRINT:  "; \
    OutputStream << std::setfill(' ') << std::setw(15) << std::left << ClassName << " "; \
    OutputStream << std::setw(26) << std::left << FunctionName  << " >  "; \
    OutputStream << std::setw(28) << std::left << VariableName << ": "; \
    OutputStream __VA_ARGS__ << NONE; \
} \


// ERROR
#define ERROR(...) \
{ \
    OutputStream << FOREGROUND_LIGHT_RED << "ERROR:  "; \
    OutputStream << std::setfill(' ') << std::setw(15) << std::left << this->GetClassName() << " "; \
    OutputStream << std::setw(26) << std::left; \
    OutputStream << std::string(__FUNCTION__) + std::string(":") + std::to_string((static_cast<long long int>(__LINE__))); \
    OutputStream << " >  "; \
    OutputStream __VA_ARGS__; \
    OutputStream << NONE << std::endl; \
} \

// WARNING
#define WARNING(...) \
{ \
    OutputStream << FOREGROUND_LIGHT_YELLOW << "WARNI:  "; \
    OutputStream << std::setfill(' ') << std::setw(15) << std::left << this->GetClassName() << " "; \
    OutputStream << std::setw(26) << std::left; \
    OutputStream << std::string(__FUNCTION__) + std::string(":") + std::to_string((static_cast<long long int>(__LINE__))); \
    OutputStream << " >  "; \
    OutputStream __VA_ARGS__; \
    OutputStream << NONE << std::endl; \
} \

// Catch Segmentation Faults
#define HERE \
{ \
    OutputStream << BACKGROUND_CYAN << "*HERE:  "; \
    OutputStream << std::setfill(' ') << std::setw(15) << std::left << this->GetClassName() << " "; \
    OutputStream << std::setw(26) << std::left << __FUNCTION__ << " >  "; \
    OutputStream << "at line: " << __LINE__; \
    OutputStream << NONE << std::endl; \
} \

// ================
// Class BaseFilter
// ================

class BaseFilter
{
    public:
        static BaseFilter * New() { return new BaseFilter; }
        inline static const char * GetFilterName() { return "BaseFilter"; }

        // Accessors, Mutators
        bool GetDebugStatus() const { return this->DebugStatus; }
        void SetDebugStatus(bool Status) { this->DebugStatus = Status; }
        void DebuggerOn() { this->DebugStatus = true; }
        void DebuggerOff() { this->DebugStatus = false; }

        bool GetDisplayStatus() const { return this->DisplayStatus; }
        void SetDisplayStatus(bool Status) { this->DisplayStatus = Status; }
        void DisplayOn() { this->DisplayStatus = true; }
        void DisplayOff() { this->DisplayStatus = false; }

        bool GetProgressStatus() const { return this->ProgressStatus; }
        void SetProgressStatus(bool Status) { this->ProgressStatus = Status; }
        void ProgressOn() { this->ProgressStatus = true; }
        void ProgressOff() { this->ProgressStatus = false; }

        unsigned int GetProgressSpeed() const { return this->ProgressSpeed; }
        void SetProgressSpeed(unsigned int Speed) { this->ProgressSpeed = Speed; }

        static bool GetLogStatus() { return LogStatus; }
        static void SetLogStatus(bool Status) { LogStatus = Status; }
        static void LoggerOn() { LogStatus = true; BaseFilter::BaseFilterOutputStream = &BaseFilter::LogFile; }
        static void LoggerOff() { LogStatus = false; BaseFilter::BaseFilterOutputStream = &std::cout; }

        static void SetLogFileName(const char *);
        static char * GetLogFileName();

    protected:
        // Constructors
        BaseFilter();
        BaseFilter(const BaseFilter & rhs);
        virtual BaseFilter * Clone() { return new BaseFilter(*this); }
        virtual ~BaseFilter();

        // Operators
        BaseFilter & operator=(const BaseFilter & rhs);

        // Member Methods
        void ProgressUpdate(unsigned int CurrentStep, unsigned int TotalSteps);

        void ProgressPrint();

        void ProgressReset();

        template<class Type>
            void DisplayFunction(
                    const char *ClassName,
                    const char *FunctionName,
                    const char *VariableName,
                    Type scalar);

        template<class Type>
            void DisplayFunction(
                    const char *ClassName,
                    const char *FunctionName,
                    const char *VariableName,
                    Type *vector,
                    unsigned int length);

        template<class Type>
            void DisplayFunction(
                    const char *ClassName,
                    const char *FunctionName,
                    const char *VariableName,
                    Type **matrix,
                    unsigned int columns,
                    unsigned int rows);
   
        template<class Type>
            int CheckSorted(Type Vector,unsigned int Length);

        // Member Data
        bool DebugStatus;
        bool DisplayStatus;
        static bool LogStatus;
        bool ProgressStatus;
        double ProgressValue; // Internal
        unsigned int ProgressSpeed;
        unsigned int ProgressInterval; // Internal
        static char *LogFileName;
        static std::ofstream & LogFile; // Internal
        static std::ostream * BaseFilterOutputStream; // Internal
};

// ==================
// Template Functions
// ==================

// ================
// Display Function
// ================

// Scalar

template<class Type>
void BaseFilter::DisplayFunction(
        const char *ClassName,
        const char *FunctionName,
        const char *VariableName,
        Type scalar)
{
    PRINT();
    OutputStream << FOREGROUND_BLUE << scalar;
    OutputStream << NONE << std::endl;
}

// Vector

template<class Type>
void BaseFilter::DisplayFunction(
        const char *ClassName,
        const char *FunctionName,
        const char *VariableName,
        Type *vector,
        unsigned int length)
{
    PRINT();
    for(unsigned int i=0; i<length; i++)
    {
        OutputStream << FOREGROUND_BLUE << vector[i];
        if(i != length-1)
        {
            OutputStream << ", ";
        }
    }
    OutputStream << NONE << std::endl;
}

// Matrix

template<class Type>
void BaseFilter::DisplayFunction(
        const char *ClassName,
        const char *FunctionName,
        const char *VariableName,
        Type **matrix,
        unsigned int columns,
        unsigned int rows)
{
    PRINT(<< std::endl);
    for(unsigned int i=0; i<rows; i++)
    {
        for(unsigned int j=0; j<columns; j++)
        {
            OutputStream << FOREGROUND_BLUE << matrix[i][j];
            if(j != columns-1)
            {
                OutputStream << ", ";
            }
        }
        OutputStream << NONE << std::endl;
    }
}

// ============
// Check Sorted
// ============

template<class Type>
int BaseFilter::CheckSorted(Type Vector,unsigned int Length)
{
    for(unsigned int i=0; i<Length-1; i++)
    {
        if(Vector[i+1]<Vector[i])
        {
            return 0;
        }
    }
    return 1;
}

#endif
