cmake_minimum_required(VERSION 2.8)

PROJECT(ApplyingVtkMarchingCubes C CXX)

set(VAMPIRTRACE_ROOT $ENV{HOME}/NaokiEto/vt-install)
set(VAMPIRTRACE_LIBRARIES ${VAMPIRTRACE_ROOT}/lib/libotfaux.so
                          ${VAMPIRTRACE_ROOT}/lib/libvt.so
                          ${VAMPIRTRACE_ROOT}/lib/libvt-hyb.so)

include_directories(${VAMPIRTRACE_ROOT}/include)
#add_definitions( -DVTRACE_PTHREAD )
add_definitions( -DVTRACE )

set(VTK_DIR /work2/VTK5.10.1-install/lib/vtk-5.10)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package (Threads)

add_executable(ApplyingVtkMarchingCubes ApplyingVtkMarchingCubes.cxx)

#SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_C_COMPILER mpicc-vt)
SET(CMAKE_C_COMPILER mpicc-vt)
SET(CMAKE_CXX_COMPILER mpicxx)
#SET(CMAKE_CXX_COMPILER vtcc++)

target_link_libraries(ApplyingVtkMarchingCubes mpi)

target_link_libraries(ApplyingVtkMarchingCubes ${VAMPIRTRACE_LIBRARIES})

target_link_libraries (ApplyingVtkMarchingCubes ${CMAKE_THREAD_LIBS_INIT})

target_link_libraries (ApplyingVtkMarchingCubes ${PTHREAD_LIBS})

if(VTK_LIBRARIES)
  target_link_libraries(ApplyingVtkMarchingCubes ${VTK_LIBRARIES})
else()
  target_link_libraries(ApplyingVtkMarchingCubes vtkHybrid)
endif()
