cmake_minimum_required(VERSION 2.8)

PROJECT(ApplyingVtkMarchingCubes)

set(VAMPIRTRACE_ROOT $ENV{HOME}/Downloads/VampirTrace-5.14.4)
set(VAMPIRTRACE_LIBRARIES ${VAMPIRTRACE_ROOT}/lib/libvt-mpi.so
                          ${VAMPIRTRACE_ROOT}/lib/libotfaux.so)

include_directories(${VAMPIRTRACE_ROOT}/include)
add_definitions( -DVTRACE )

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(ApplyingVtkMarchingCubes ApplyingVtkMarchingCubes.cxx)

SET(CMAKE_C_COMPILER mpicc-vt)
SET(CMAKE_CXX_COMPILER mpicxx)

target_link_libraries(ApplyingVtkMarchingCubes ${VAMPIRTRACE_LIBRARIES})

target_link_libraries(ApplyingVtkMarchingCubes mpi)

if(VTK_LIBRARIES)
  target_link_libraries(ApplyingVtkMarchingCubes ${VTK_LIBRARIES})
else()
  target_link_libraries(ApplyingVtkMarchingCubes vtkHybrid)
endif()
