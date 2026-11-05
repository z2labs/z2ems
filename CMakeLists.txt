
## usage: cmake -DCMAKE_INSTALL_PREFIX=<install/path> .
## use -DBUILD_APPCSXCAD=NO to not build QCSXCAD and AppCSXCAD

# In CMake 4, 3.10 is deprecated and 3.5 has been removed.
# use 3.0...3.10 so all of these versions are acceptable as min. version.
# https://cmake.org/cmake/help/latest/command/cmake_minimum_required.html
cmake_minimum_required(VERSION 3.0...3.10)

PROJECT(openEMS-Project CXX)

if (CMAKE_TOOLCHAIN_FILE)
    message(STATUS "Using CMAKE_TOOLCHAIN_FILE: " ${CMAKE_TOOLCHAIN_FILE})
endif()

if (CMAKE_INSTALL_PREFIX)
    message(STATUS "Using CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
endif()

if (NOT DEFINED BUILD_APPCSXCAD)
    set(BUILD_APPCSXCAD "YES")
endif()

message(STATUS "Build AppCSXCAD: " ${BUILD_APPCSXCAD})

IF(EXISTS ${PROJECT_SOURCE_DIR}/localConfig.cmake)
   message(STATUS "include local config" )
   include(${PROJECT_SOURCE_DIR}/localConfig.cmake)
ENDIF()

include(ExternalProject)

# MSVC: binaries built with the 14.40+ toolset (VS 2022 17.10 and later) need an equally new msvcp140.dll on the
# target machine, because std::mutex got a constexpr constructor there; on a machine with an older VC++ runtime
# the first mutex lock faults (0xC0000005 in MSVCP140.dll / 0xC0000409 fail-fast). In openEMS that is the first
# barrier of the multithreaded engine, so the single-threaded engines still run and the crash looks like a CPU
# problem. _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR restores the runtime-initialised mutex, which works with every
# msvcp140.dll. (Shipping the toolset's own CRT next to the executables is the other half: scripts/win_bundle.ps1.)
set(MSVC_COMPAT_ARGS "")
if (MSVC)
    set(MSVC_COMPAT_ARGS "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} /D_DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR")
endif()

# build fparser
ExternalProject_Add( fparser
  SOURCE_DIR  ${PROJECT_SOURCE_DIR}/fparser
  CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${MSVC_COMPAT_ARGS}
)

# build CSXCAD
ExternalProject_Add( CSXCAD 
  DEPENDS     fparser
  SOURCE_DIR  ${PROJECT_SOURCE_DIR}/CSXCAD
  CMAKE_ARGS  -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DFPARSER_ROOT_DIR=${CMAKE_INSTALL_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${MSVC_COMPAT_ARGS}
)

# CUDA engine of openEMS (--engine=cuda): -DWITH_CUDA=ON, optionally -DCMAKE_CUDA_ARCHITECTURES="61;75;86;89"
# (61 = Pascal, GTX 10xx: needs a CUDA 12.x toolkit, CUDA 13 dropped Pascal)
option(WITH_CUDA "build the CUDA engine of openEMS" OFF)
set(OPENEMS_CUDA_ARGS "")
if (WITH_CUDA)
    set(OPENEMS_CUDA_ARGS -DWITH_CUDA=ON)
    if (DEFINED CMAKE_CUDA_ARCHITECTURES)
        string(REPLACE ";" "|" _cuda_archs "${CMAKE_CUDA_ARCHITECTURES}")
        list(APPEND OPENEMS_CUDA_ARGS "-DCMAKE_CUDA_ARCHITECTURES=${_cuda_archs}")
    endif()
endif()
message(STATUS "Build the CUDA engine: " ${WITH_CUDA})

# build openEMS
ExternalProject_Add( openEMS
  DEPENDS     fparser CSXCAD
  SOURCE_DIR  ${PROJECT_SOURCE_DIR}/openEMS
  LIST_SEPARATOR |
  CMAKE_ARGS  -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DFPARSER_ROOT_DIR=${CMAKE_INSTALL_PREFIX} -DCSXCAD_ROOT_DIR=${CMAKE_INSTALL_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${MSVC_COMPAT_ARGS} ${OPENEMS_CUDA_ARGS}
)

if (${BUILD_APPCSXCAD})
    # build QCSXCAD
    ExternalProject_Add( QCSXCAD
    DEPENDS     CSXCAD
    SOURCE_DIR  ${PROJECT_SOURCE_DIR}/QCSXCAD
    CMAKE_ARGS  -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DCSXCAD_ROOT_DIR=${CMAKE_INSTALL_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${MSVC_COMPAT_ARGS}
    )

    # build AppCSXCAD
    ExternalProject_Add( AppCSXCAD
    DEPENDS     QCSXCAD
    SOURCE_DIR  ${PROJECT_SOURCE_DIR}/AppCSXCAD
    CMAKE_ARGS  -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DCSXCAD_ROOT_DIR=${CMAKE_INSTALL_PREFIX} -DQCSXCAD_ROOT_DIR=${CMAKE_INSTALL_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${MSVC_COMPAT_ARGS}
    )
endif()
