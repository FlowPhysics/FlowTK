#
# =====================================================================================
#
#       Filename:  CMakeLists.txt
#
#    Description:  Cmake file for Pipeline
#
#        Version:  1.0
#        Created:  10/25/2012 21:00:40 PM
#       Revision:  none
#       Compiler:  gcc
#
#         Author:  Siavash Ameli
#   Organization:  University of California, Berkeley
#
# =====================================================================================

# ===================
# Directory Structure
#
# Project
#   ^
#   |---> bin
#   |---> build
#   |---> config
#   |---> doc
#   |---> install
#   |---> include
#   |---> src
#   |---> test
#   |---> lib
#   |---> plugin
#
# ===================

cmake_minimum_required(VERSION 2.8)
project(flowVP CXX)

# ======================
# VTK/ParaView Inclusion
# ======================

# ParaView plugin option
set(BUILD_PARAVIEW_PLUGIN ON CACHE BOOL "Build ParaView plugin?")
# set(BUILD_PARAVIEW_PLUGIN OFF CACHE BOOL "Build ParaView plugin?")

# Avoid confliction of Paraview and VTK inclusion
if(BUILD_PARAVIEW_PLUGIN)
    # include ParaView instead of VTK
    find_package(ParaView REQUIRED)
    if(ParaView_FOUND)
        include(${PARAVIEW_USE_FILE})
    else(ParaView_FOUND)
        message(FATAL_ERROR "ParaView not found.")
    endif(ParaView_FOUND)
else(BUILD_PARAVIEW_PLUGIN)
    # include VTK
    find_package(VTK REQUIRED)
    if(VTK_FOUND)
        include(${VTK_USE_FILE})
    else(VTK_FOUND)
        message(FATAL_ERROR "VTK not found.")
    endif(VTK_FOUND)
endif(BUILD_PARAVIEW_PLUGIN)

# ==============
# Compiler Flags
# ==============

# Initialize CXX FLAGS
set(CMAKE_CXX_FLAGS                "-Wall -std=c++0x ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# Compiler-specific C++11 activation
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.4 OR GCC_VERSION VERSION_EQUAL 4.4))
        message(FATAL_ERROR "${PROJECT_NAME} rRequires g++ 4.4 or greater.")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif()

# ===================
# Directory Structure
# ===================

# Headers
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${PROJECT_INCLUDE_DIR})
set(TEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
include_directories(${TEST_INCLUDE_DIR})

# Library
set(PROJECT_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
aux_source_directory(${PROJECT_LIBRARY_DIR} LIBRARY_SOURCES)
add_library(FILTERS_LIBRARY STATIC ${LIBRARY_SOURCES})

# Source
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
aux_source_directory(${PROJECT_SOURCE_DIR} SOURCES)
set(MAIN_EXECUTABLE_FILENAME "exe")
add_executable(${MAIN_EXECUTABLE_FILENAME} WIN32 MACOSX_BUNDLE ${SOURCES})
if(VTK_LIBRARIES)
    target_link_libraries(${MAIN_EXECUTABLE_FILENAME} FILTERS_LIBRARY ${VTK_LIBRARIES})
else(VTK_LIBRARIES)
    target_link_libraries(${MAIN_EXECUTABLE_FILENAME} FILTERS_LIBRARY vtkHybrid)
endif(VTK_LIBRARIES)

# Test
set(PROJECT_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
aux_source_directory(${PROJECT_TEST_DIR} TESTS)
foreach (FILE ${TESTS})
    get_filename_component(TEST_EXECUTABLE_NAME ${FILE} NAME_WE)
    add_executable(${TEST_EXECUTABLE_NAME} WIN32 MACOSX_BUNDLE ${FILE})
    if(VTK_LIBRARIES)
        target_link_libraries(${TEST_EXECUTABLE_NAME} FILTERS_LIBRARY ${VTK_LIBRARIES})
    else(VTK_LIBRARIES)
        target_link_libraries(${TEST_EXECUTABLE_NAME} FILTERS_LIBRARY vtkHybrid)
    endif(VTK_LIBRARIES)
endforeach(FILE)

# ================
# ParaView Plugins
# ================

# TODO: add documetation for plugins

if(BUILD_PARAVIEW_PLUGIN)
    set(PROJECT_PLUGIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/plugin)
    set(PROJECT_PLUGIN_DOCUMENTATION_DIR ${CMKAE_CURRENT_SOURCE_DIR}/doc)
    file(GLOB PLUGIN_XML_FILES ${PROJECT_PLUGIN_DIR}/*[!GUI].xml)
    set(PluginVersion "1.0")

    # Loop over each plugin
    foreach (XML_FILE ${PLUGIN_XML_FILES})
        get_filename_component(FilenameWithoutExtension ${XML_FILE} NAME_WE)
        set(CXX_FILE "${PROJECT_LIBRARY_DIR}/${FilenameWithoutExtension}.cxx")
        set(GUI_FILE "${PROJECT_PLUGIN_DIR}/${FilenameWithoutExtension}GUI.xml")
        set(ICN_FILE "${PROJECT_PLUGIN_DIR}/${FilenameWithoutExtension}ICN.qrc")
        set(DOC_FILE "${PROJECT_PLUGIN_DIR}/${FilenameWithputExtention}.html")

        if(EXISTS ${GUI_FILE} AND EXISTS ${ICN_FILE} AND EXISTS)
            # Plugin with GUI, Icon and Doc
            Add_PARAVIEW_PLUGIN(
                ${FilenameWithoutExtension}
                ${PluginVersion}
                SERVER_MANAGER_XML      ${XML_FILE}
                DOCUMENTATION_DIR       ${PROJECT_PLUGIN_DOCUMENTATION_DIR}
                GUI_RESOURCE_FILES      ${GUI_FILE}
                GUI_RESOURCES           ${ICN_FILE}
                SERVER_MANAGER_SOURCES  ${CXX_FILE})

        elseif(EXISTS ${GUI_FILE} AND EXISTS ${ICN_FILE})
            # Plugin with GUI and Icon
            Add_PARAVIEW_PLUGIN(
                ${FilenameWithoutExtension}
                ${PluginVersion}
                SERVER_MANAGER_XML      ${XML_FILE}
                DOCUMENTATION_DIR       ${PROJECT_PLUGIN_DOCUMENTATION_DIR}
                GUI_RESOURCE_FILES      ${GUI_FILE}
                GUI_RESOURCES           ${ICN_FILE}
                SERVER_MANAGER_SOURCES  ${CXX_FILE})

        elseif(EXISTS ${GUI_FILE})
            # Plugin with GUI file
            ADD_PARAVIEW_PLUGIN(
                ${FilenameWithoutExtension}
                ${PluginVersion}
                SERVER_MANAGER_XML      ${XML_FILE}
                GUI_RESOURCE_FILES      ${GUI_FILE}
                SERVER_MANAGER_SOURCES  ${CXX_FILE})

        else()
            # Plugin without GUI file
            ADD_PARAVIEW_PLUGIN(
                ${FilenameWithoutExtension}
                ${PluginVersion}
                SERVER_MANAGER_XML      ${XML_FILE}
                SERVER_MANAGER_SOURCES  ${CXX_FILE})
        endif()

    endforeach(XML_FILE)
endif(BUILD_PARAVIEW_PLUGIN)

# ==================
# Output Directories
# ==================

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Directory for all executables.")
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Directory for all libraries.")

# ======================
# Generate Documentation
# ======================

# TODO: doxygen
set(BUILD_DOC ON CACHE BOOL "Build documentation?")

if(BUILD_DOC)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        set(DOXYGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doc)
        set(DOXYGEN_FILE ${DOXYGEN_DIR}/Doxyfile.in)
        # Some stuff here
    else(DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen not found.")
    endif(DOXYGEN_FOUND)
endif(BUILD_DOC)

# =========
# Build GUI
# =========

# TODO: Qt4
set(BUILD_GUI ON CACHE BOOL "Build GUI?")

find_package(Qt4 REQUIRED)
if(QT4_FOUND)
    include(${QT_USE_FILE})
    # refer to page 49-50 of cmake.
else(QT4_FOUND)
    message(FATAL_ERROR "Qt4 not found.")
endif(QT4_FOUND)

# ================
# Enable Debugging
# ================

# TODO: Still it is not shown in GUI.
if(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Options are: empty, Debug,Release, RelWithDebInfo,MinSizeRel")
else(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Options are: empty, Debug,Release, RelWithDebInfo,MinSizeRel" FORCE)
endif(DEFINED CMAKE_BUILD_TYPE)

# ======================
# Windows Considerations
# ======================

if(WIN32)
    set(CMAKE_EXECUTABLE_SUFFIX ".exe")
endif(WIN32)

# ========================
# Get Number Of Processors
# ========================

# Determine number of processors on host system
if(NOT DEFINED NUMBER_OF_HOST_PROCESSORS)
    set(NUMBER_OF_HOST_PROCESSORS 0)

    # Linux
    if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        set(cpuinfo_file "/proc/cpuinfo")
        if(EXISTS "${cpuinfo_file}")
            file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
            list(LENGTH procs NUMBER_OF_HOST_PROCESSORS)
        endif()

    # Windows:
    elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
        set(NUMBER_OF_HOST_PROCESSORS "$ENV{NUMBER_OF_PROCESSORS}")

    # Apple
    elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
        find_program(cmd_sys_pro "system_profiler")
        if(cmd_sys_pro)
            execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
            string(REGEX REPLACE "^.*Total Number Of Cores: ([0-9]+).*$" "\\1" NUMBER_OF_HOST_PROCESSORS "${info}")
        endif(cmd_sys_pro)
    endif()
endif()

# Determine number of processors on target system
if(NOT DEFINED NUMBER_OF_TARGET_PROCESSORS)
    set(NUMBER_OF_TARGET_PROCESSORS 0)

    # Linux
    if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        set(cpuinfo_file "/proc/cpuinfo")
        if(EXISTS "${cpuinfo_file}")
            file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
            list(LENGTH procs NUMBER_OF_TARGET_PROCESSORS)
        endif()

    # Windows:
    elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
        set(NUMBER_OF_TARGET_PROCESSORS "$ENV{NUMBER_OF_PROCESSORS}")

    # Apple
    elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
        find_program(cmd_sys_pro "system_profiler")
        if(cmd_sys_pro)
            execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
            string(REGEX REPLACE "^.*Total Number Of Cores: ([0-9]+).*$" "\\1" NUMBER_OF_TARGTE_PROCESSORS "${info}")
        endif(cmd_sys_pro)
    endif()
endif()

# ======
# CCache
# ======

# TODO: make this for Apple and windows

# Avoid using ccache if Qt is enabled.
if(NOT BUILD_PARAVIEW_PLUGIN AND NOT BUILD_GUI)
    find_program(CCACHE "ccache" PATHS "/usr/bin")
    if(NOT CCACHE MATCHES CCACHE-NOTFOUND)
        set(CMAKE_C_COMPILER ${CCACHE})
        set(CMAKE_CXX_COMPILER ${CCACHE})
        set(CMAKE_C_FLAGS "gcc ${CMAKE_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "g++ ${CMAKE_CXX_FLAGS}")
    else(NOT CCACHE MATCHES CCACHE-NOTFOUND)
        message("ccache not found.")
    endif(NOT CCACHE MATCHES CCACHE-NOTFOUND)
endif(NOT BUILD_PARAVIEW_PLUGIN AND NOT BUILD_GUI)

# =====================
# Build Parallelization
# =====================

if(${NUMBER_OF_HOST_PROCESSORS} GREATER 1)
    math(EXPR NUMBER_OF_EFFECTIVE_HOST_JOBS "${NUMBER_OF_HOST_PROCESSORS} + 1")
    if(${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
        set(CMAKE_MAKE_PROGRAM "${CMAKE_MAKE_PROGRAM} -j ${NUMBER_OF_EFFECTIVE_HOST_JOBS}")
    endif(${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
endif(${NUMBER_OF_HOST_PROCESSORS} GREATER 1)

# ===============
# Cross Compiling
# ===============

# TODO; find a way to put -DCMAKE_TOOLCHAIN_FILE=../Toolchain-windows.cmake
#       into CMakeLists.txt fie.

# Check mingw32 is installed
find_program(MINGW NAME "i586-mingw32msvc-g++" PATHS "/usr/bin")
if(NOT MINGW MATCHES MINGW-NOTFOUND)
    # add cross compilation here
else()
    message(FATAL_ERROR "mingw32 not found")
endif()

# ============
# Installation
# ============

# TODO: add support for "make uninstall", refer to Cmake FAQ in wiki.

option(INSTALL_PROJECT "Install project on host system." ON)
if(${INSTALL_PROJECT})
    set(INSTALL_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/install CACHE PATH "Path to install project.")
    install(TARGETS ${MAIN_EXECUTABLE_FILENAME} DESTINATION ${INSTALL_PROJECT_PATH})
endif(${INSTALL_PROJECT})

# Note: Cpack also uses the same install path.

# ================
# Generate Package
# ================

option(GENERATE_PACKAGE "Generate system dependent package" ON)

if(${GENERATE_PACKAGE})
    include(InstallRequiredSystemLibraries)

    # Custom Variables
    set(AUTHOR "Siavash Ameli")
    set(PACKAGE_NAME "flowVP")

    # Cpack variables
    set(CPACK_PACKAGE_NAME ${PACKAGE_NAME})
    set(CPACK_PACKAGE_CONTACT "sameli@iit.edu")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
    set(CPACK_PACKAGE_VERSION_MAJOR 1)
    set(CPACK_PACKAGE_VERSION_MINOR 0)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "flow Visualization Pipeline/Plugin")
    set(CPACK_PACKAGE_FILE_NAME ${PACKAGE_NAME})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")

    # Package install path
    if(NOT DEFINED INSTALL_PROJECT_PATH)
        set(CPACK_SET_DESTDIR true)
        set(CPACK_INSTALL_PREFIX ${CMAKE_CURRENT_PROJECT}/install)
    endif(NOT DEFINED INSTALL_PROJECT_PATH)

    option(BINARY_PACKAGE "Generate binary package" ON)

    # Operating System
    # 1- Linux
    if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        if(${BINARY_PACKAGE})
            set(CPACK_GENERATOR "DEB")
            set(CPACK_DEBIAN_PACKAGE_MAINTAINER ${AUTHOR})
            set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
        else(${BINARY_PACKAGE})
            set(CPACK_GENERTOR "SOURCE_TGZ")
        endif(${BINARY_PACKAGE})

    # 2- Windows
    elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
        set(CPACK_NSIS_CONTACT ${AUTHOR})
        set(CPACK_NSIS_PACKAGE_NAME ${PACKAGE_NAME})
        if(${BINARY_PACKAGE})
            set(CPACK_GENERATOR "NSIS")
        else(${BINARY_PACKAGE})
            set(CPACK_GENERATOR "SOURCE_ZIP")
        endif(${BINAY_PACKAGE})

    # 3- Apple
    elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
        set(CPACK_GENERATOR "Drag and Drop")    
    endif()

    include(CPack)

endif(${GENERATE_PACKAGE})

# =====
# CTest
# =====

option(TEST_PROJECT "Test the project using CTest" ON)

if(${TEST_PROJECT})
    include(CTest)
    ENABLE_TESTING()
    ADD_TEST(test exe)

    # Parallelizatoin
    if(NUMBER_OF_TARGET_PROCESSORS)
        set(CTEST_BUILD_FLAGS "-j ${NUMBER_OF_TARGET_PROCESSORS}")
    endif(NUMBER_OF_TARGET_PROCESSORS)
endif(${TEST_PROJECT})

# =================
# Print Information
# =================

option(PRINT_VARIABLES "Print useful information about make process." OFF)

if(${PRINT_VARIABLES})
    message("CMAKE_BINARY_DIR:           ${CMAKE_BINARY_DIR}")
    message("CMAKE_CURRENT_BINARY_DIR:   ${CMAKE_CURRENT_BINARY_DIR}")
    message("CMAKE_CURRENT_LIST_FILE:    ${CMAKE_CURRENT_LIST_FILE}")
    message("CMAKE_CURRENT_LIST_DIR:     ${CMAKE_CURRENT_LIST_DIR}")
    message("CMAKE_CURRENT_SOURCE_DIR:   ${CMAKE_CURRENT_SOURCE_DIR}")
    message("CMAKE_SOURCE_DIR:           ${CMAKE_SOURCE_DIR}")
    message("PROJECT_BINARY_DIR:         ${PROJECT_BINARY_DIR}")
    message("PROJECT_SOURCE_DIR:         ${PROJECT_SOURCE_DIR}")
    message("EXECUTABLE_OUTPUT_PATH:     ${EXECUTABLE_OUTPUT_PATH}")
    message("LIBRARY_OUTPUT_PATH:        ${LIBRARY_OUTPUT_PATH}")
endif(${PRINT_VARIABLES})

# =======================
# Hide Variables from GUI
# =======================

mark_as_advanced(FORCE CMAKE_INSTALL_PREFIX VTK_DIR CCACHE)
