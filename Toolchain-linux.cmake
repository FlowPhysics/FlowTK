# =====================================================================================
#
#       Filename:  Toolchain-linux.cmake
#
#    Description:  Toolchain file for cross compiling to embedded linux system target
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

# =============================
# Host system:   Windows
# Target system: Embedded Linux
# =============================

# Target system name
set(CMAKE_SYSTEM_NAME Linux)

# Compilers for C and C++
set(CMAKE_C_COMPILER   /usr/bin/gcc)
set(CMAKE_CXX_COMPILER /usr/bin/g++)

# Target environment
set(CMAKE_FIND_ROOT_PATH /usr)

# set defaults of find_xyz()
set(CMKAE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
