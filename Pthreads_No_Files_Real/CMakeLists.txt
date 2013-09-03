
cmake_minimum_required(VERSION 2.8)

PROJECT(ApplyingVtkMarchingCubes C CXX)

set(VAMPIRTRACE_ROOT /work2/vt-install)
set(VAMPIRTRACE_LIBRARIES ${VAMPIRTRACE_ROOT}/lib/libotfaux.so
                          ${VAMPIRTRACE_ROOT}/lib/libvt.so
                          ${VAMPIRTRACE_ROOT}/lib/libvt-mt.so)

include_directories(${VAMPIRTRACE_ROOT}/include)
add_definitions( -DVTRACE_PTHREAD )

set(VTK_DIR /work2/VTK5.10.1-install/lib/vtk-5.10)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package (Threads)

add_executable(ApplyingVtkMarchingCubes ApplyingVtkMarchingCubes.cxx)

target_link_libraries(ApplyingVtkMarchingCubes ${VAMPIRTRACE_LIBRARIES})

target_link_libraries (ApplyingVtkMarchingCubes ${CMAKE_THREAD_LIBS_INIT})

target_link_libraries (ApplyingVtkMarchingCubes ${PTHREAD_LIBS})

if(VTK_LIBRARIES)
  target_link_libraries(ApplyingVtkMarchingCubes ${VTK_LIBRARIES})
else()
  target_link_libraries(ApplyingVtkMarchingCubes vtkHybrid)
endif()
