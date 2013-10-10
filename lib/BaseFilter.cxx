/*
 * =====================================================================================
 *
 *       Filename:  BaseFilter.cxx
 *
 *    Description:  Base Filter
 *
 *        Version:  1.0
 *        Created:  09/12/2012 03:07:40 AM
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

#include <BaseFilter.h>
#include <cstring>    // for strcpy

// ======
// Macros
// ======

#define MAX_CHAR_LENGTH 512

// ================
// Static variables
// ================

static std::ofstream TemporaryLogFile(NULL);
std::ofstream & BaseFilter::LogFile = TemporaryLogFile;

// Default Log status
//std::ostream * BaseFilter::BaseFilterOutputStream = & BaseFilter::LogFile;  // Log on
std::ostream * BaseFilter::BaseFilterOutputStream = & std::cout;          // Log off

char * BaseFilter::LogFileName = new char[512];
bool BaseFilter::LogStatus = false;

// ===========
// Constructor
// ===========

BaseFilter::BaseFilter()
{
    // Member Data
    this->DebuggerOff();
    this->DisplayOff();
    this->LoggerOff();
    this->ProgressOff();
    this->SetProgressSpeed(1);
    this->ProgressMessage = NULL;

    // Internal Member Data
    this->ProgressValue = 0;
    this->ProgressInterval = 0;
}

// ===============
// Copy Custructor
// ===============

BaseFilter::BaseFilter(const BaseFilter & rhs)
{
    // Copy Member data
    this->DebugStatus = rhs.GetDebugStatus();
    this->DisplayStatus = rhs.GetDisplayStatus();
    this->LogStatus = rhs.GetLogStatus();
    this->ProgressStatus = rhs.GetProgressStatus();
    this->ProgressSpeed = rhs.GetProgressSpeed();
    this->SetProgressMessage(rhs.GetProgressMessage());
}

// ==========
// Destructor
// ==========

BaseFilter::~BaseFilter()
{
    // Log File
    LogFile.close();
    delete [] BaseFilter::LogFileName;
    BaseFilter::LogFileName = NULL;
    BaseFilter::BaseFilterOutputStream = NULL;

    // Progress Message
    if(this->ProgressMessage != NULL)
    {
        delete [] this->ProgressMessage;
        this->ProgressMessage = NULL;
    }
}

// ========
// Operator
// ========

BaseFilter & BaseFilter::operator=(const BaseFilter & rhs)
{
    if(this == &rhs)
    {
        return *this;
    }
    this->DebugStatus = rhs.GetDebugStatus();
    this->DisplayStatus = rhs.GetDisplayStatus();
    this->LogStatus = rhs.GetLogStatus();
    this->ProgressStatus = rhs.GetProgressStatus();
    this->ProgressSpeed = rhs.GetProgressSpeed();
    return *this;
}

// ====================
// Get Progress Message
// ====================

inline const char * BaseFilter::GetProgressMessage() const
{
    if(this->ProgressMessage == NULL)
    {
        return "";
    }
    else
    {
        return this->ProgressMessage;
    }
}

// ====================
// Set Progress Message
// ====================

inline void BaseFilter::SetProgressMessage(const char *Message)
{
    // Delete previous message
    if(this->ProgressMessage != NULL)
    {
        delete [] this->ProgressMessage;
    }

    // Allocate memory
    this->ProgressMessage = new char[MAX_CHAR_LENGTH];

    // Set Message
    std::strcpy(this->ProgressMessage,Message);
}

// =================
// Set Log File Name
// =================

void BaseFilter::SetLogFileName(const char * FileName)
{
    std::strcpy(BaseFilter::LogFileName,FileName);
    LogFile.open(BaseFilter::LogFileName,std::ios_base::out);
    if(!LogFile.is_open())
    {
        std::cerr << FOREGROUND_RED << "ERROR: BaseFilter >> SetLogFileName >> LogFile is not open." << NONE << std::endl;
    }
}

// =================
// Get Log File Name
// =================

char * BaseFilter::GetLogFileName()
{
    return BaseFilter::LogFileName;
}

// ===============
// Progress Update
// ===============

void BaseFilter::ProgressUpdate(unsigned int CurrentStep, unsigned int TotalSteps)
{
    // Initialize Progress Interval
    if(this->ProgressInterval == 0)
    {
        this->ProgressInterval = floor(TotalSteps/(100*this->ProgressSpeed));
        this->ProgressInterval = (this->ProgressInterval >= 1 ? this->ProgressInterval : 1);
    }

    // Update Progress value
    if(CurrentStep % this->ProgressInterval == 0 || CurrentStep == TotalSteps-1)
    {
        this->ProgressValue = 100*static_cast<double>(CurrentStep+1)/TotalSteps;
    }

    // Print Progress
    if(this->GetProgressStatus() == true)
    {
        this->ProgressPrint();
    }
}

// ==============
// Progress Print
// ==============

void BaseFilter::ProgressPrint()
{
    std::cout.setf(ios::fixed,ios::floatfield);
    std::cout.precision(2);
    std::cout << "\r" << FOREGROUND_YELLOW << "PROGR:  ";
    std::cout << std::setfill(' ') << std::setw(OutputStreamSecondColumnWidth) << std::left << this->GetFilterName() << " ";
    std::cout << std::setw(OutputStreamThirdColumnWidth) << std::left << this->GetProgressMessage();  
    std::cout << " >  Progress: ";
    std::cout << std::setprecision(2) << std::right << std::setw(6) << this->ProgressValue << " %";
    std::cout << NONE;
    if(this->ProgressValue == 100)
    {
        std::cout << FOREGROUND_YELLOW << "  >  Done!";
        std::cout << NONE << std::endl;
    }
}

// ==============
// Progress Reset
// ==============

void BaseFilter::ProgressReset()
{
    this->ProgressValue = 0;
    this->ProgressInterval = 0;
    
    // Progress Message
    if(this->ProgressMessage != NULL)
    {
        delete [] this->ProgressMessage;
        this->ProgressMessage = NULL;
    }
}
