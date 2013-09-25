# =====================================================================================
#
#       Filename:  Toolchain-windows.cmake
#
#    Description:  Toolchain file for cross compiling to windows system target
#
#        Version:  1.0
#        Created:  12/11/2012 15:01:32 PM
#       Revision:  none
#       Compiler:  gcc
#
#         Author:  Siavash Ameli
#   Organization:  University of California, Berkeley
#
# =====================================================================================

# ======================
# Host system:   Linux
# Target system: Windows
# ======================

# Target system name
set(CMAKE_SYSTEM_NAME Windows)

# Compilers for C and C++
set(CMAKE_C_COMPILER   i586-mingw32msvc-gcc)
set(CMAKE_CXX_COMPILER i586-mingw32msvc-g++)
set(CMAKE_RC_COMPILER  i586-mingw32msvc-windres)

# Target environment
set(CMAKE_FIND_ROOT_PATH /usr/i586-mingw32msvc /home/sia/programs/vtk/vtk-5.10.0-build)

# set defaults of find_xyz()
set(CMKAE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
