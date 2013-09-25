cmake_minimum_required(VERSION 2.8)
project(TestIndividualFilters CXX)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# C++11
# Initialize CXXFLAGS.
set(CMAKE_CXX_FLAGS                "-Wall -std=c++0x ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.4 OR GCC_VERSION VERSION_EQUAL 4.4))
        message(FATAL_ERROR "${PROJECT_NAME} Requires g++ 4.4 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

# Set paths
set(PROJECT_DIR "/home/sia/code/projvtk/Filter/Pipeline")
set(HEADERS_DIR ${PROJECT_DIR}/include)
set(LIBRARY_DIR ${PROJECT_DIR}/lib)

# Include Headers
include_directories(${HEADERS_DIR})
include_directories(${LIBRARY_DIR})

# Executable Name
# set(EXECUTABLE_NAME "Reader")
set(EXECUTABLE_NAME "Seed")
# set(EXECUTABLE_NAME "Cache")
# set(EXECUTABLE_NAME "Interpolator")
# set(EXECUTABLE_NAME "FlowMap")
# set(EXECUTABLE_NAME "Deformation")
# set(EXECUTABLE_NAME "LCS")
# set(EXECUTABLE_NAME "Consumer")

# Filter Filename
set(FILTER_FILENAME ${LIBRARY_DIR}/${EXECUTABLE_NAME}.cxx)

# Base Filter Filename
set(BASE_FILTER_FILENAME ${LIBRARY_DIR}/BaseFilter.cxx)

# Add executables
# add_executable(${EXECUTABLE_NAME} ${FILTER_FILENAME} ${BASE_FILTER_FILENAME})

# Add Library
add_library(FILTER_LIBRARY ${FILTER_FILENAME} ${BASE_FILTER_FILENAME})

# Target libraries
# if(VTK_LIBRARIES)
#     target_link_libraries(${EXECUTABLE_NAME} ${VTK_LIBRARIES})
# else(VTK_LIBRARIES)
#     target_link_libraries(${EXECUTABLE_NAME} vtkHybrid)
# endif(VTK_LIBRARIES)

message(" ")
message("=============================================")
message("THIS IS INDIVIDUAL TEST BUILD FOR THE PROJECT")
message("=============================================")
message(" ")
