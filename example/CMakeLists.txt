cmake_minimum_required(VERSION 3.13...3.20)
project(p4estExamples LANGUAGES C)

include(CTest)
include(CheckIncludeFile)
include(CheckSymbolExists)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules)

option(mpi "use MPI" on)

# --- find external libraries
if(mpi)
  find_package(MPI COMPONENTS C REQUIRED)
  if(NOT MPIEXEC_EXECUTABLE)
    message(FATAL_ERROR "MPIEXEC not found")
  endif()
endif(mpi)

find_package(ZLIB REQUIRED)

# --- find our library
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

# SC could have been installed with p4est,
# or SC could be installed by itself
find_package(SC HINTS ${P4EST_ROOT} PATHS ${PROJECT_SOURCE_DIR}/../build REQUIRED)

find_package(P4EST COMPONENTS P8EST PATHS ${PROJECT_SOURCE_DIR}/../build)
if(NOT P4EST_FOUND)
  find_package(P4EST REQUIRED PATHS ${PROJECT_SOURCE_DIR}/../build)
endif()

# --- get system capabilities
check_include_file(arpa/inet.h P4EST_HAVE_ARPA_INET_H)
check_include_file(netinet/in.h P4EST_HAVE_NETINET_IN_H)
if(WIN32 AND NOT P4EST_HAVE_ARPA_INET_H AND NOT P4EST_HAVE_NETINET_IN_H)
  check_include_file(Winsock2.h P4EST_HAVE_WINSOCK2_H)
  set(WINSOCK_LIBRARIES wsock32 ws2_32) # Iphlpapi
endif()

if(NOT (P4EST_HAVE_ARPA_INET_H OR P4EST_HAVE_NETINET_IN_H OR P4EST_HAVE_WINSOCK2_H))
  message(FATAL_ERROR "A networking library was not found.")
endif()

check_include_file(unistd.h P4EST_HAVE_UNISTD_H)
if(P4EST_HAVE_UNISTD_H)
  check_include_file(getopt.h P4EST_HAVE_GETOPT_H)
endif()


cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
if(Ncpu LESS 2)
  include(ProcessorCount)
  ProcessorCount(n)
  if(n GREATER Ncpu)
    set(Ncpu ${n})
  endif()
endif()

# --- helper functions

function(p4est_example name files arg1 arg2)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P4EST SC::SC)
if(MPI_FOUND)
  target_link_libraries(${name} PRIVATE MPI::MPI_C)
  add_test(NAME p4est:${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:${name}> ${arg1} ${arg2})
else()
  add_test(NAME p4est:${name} COMMAND $<TARGET_FILE:${name}> ${arg1} ${arg2})
endif()
set_tests_properties(p4est:${name} PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi)

endfunction(p4est_example)


function(p8est_example name files arg1 arg2)

add_executable(${name} ${files})
target_link_libraries(${name} PRIVATE P4EST::P8EST SC::SC)
if(MPI_FOUND)
  target_link_libraries(${name} PRIVATE MPI::MPI_C)
  add_test(NAME p8est:${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:${name}> ${arg1} ${arg2})
else()
  add_test(NAME p8est:${name} COMMAND $<TARGET_FILE:${name}> ${arg1} ${arg2})
endif()
set_tests_properties(p8est:${name} PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi)

endfunction(p8est_example)

# --- setup examples

# p4est_example(mesh2 mesh/mesh2.c unit 4)  # mesh/mesh2.c:291:9: error: too many arguments to function ‘p4est_mesh_quadrant_cumulative’
# p4est_example(points2 points/points2.c unit 4)  # missing .pts input files

p4est_example(simple2 simple/simple2.c unit 4)
if(TARGET P4EST::P8EST)
  p8est_example(simple3 simple/simple3.c unit 4)
endif()


if(P4EST_HAVE_GETOPT_H)

foreach(n bricks timings)
  p4est_example(${n}2 timings/${n}2.c "" "")
  if(TARGET P4EST::P8EST)
    p8est_example(${n}3 timings/${n}3.c "" "")
  endif()
endforeach()

# optparse not recognizing -l --level
# p4est_example(loadconn2 timings/loadconn2.c "-l 4" "")

set(n spheres)
p4est_example(${n}2 "${n}/${n}2.c;${n}/p4est_${n}.c" "" "")
if(TARGET P4EST::P8EST)
  p8est_example(${n}3 "${n}/${n}3.c;${n}/p8est_${n}.c" "" "")
endif()

set(n particles)
p4est_example(${n}2 ${n}/${n}2.c "" "")
if(TARGET P4EST::P8EST)
  p8est_example(${n}3 ${n}/${n}3.c "" "")
endif()

endif(P4EST_HAVE_GETOPT_H)

foreach(i 1 3 4)
  set(n p4est_step${i})
  p4est_example(${n} steps/${n}.c "" "")

  if(TARGET P4EST::P8EST)
    set(n p8est_step${i})
    p8est_example(${n} steps/${n}.c "" "")
  endif()
endforeach()

# optparse not recognizing -l --level
# p4est_example(p4est_step2 steps/p4est_step2.c ${CMAKE_CURRENT_SOURCE_DIR}/steps/hole_2d_cubit.inp 4)
# set_tests_properties(p4est:p4est_step2 PROPERTIES TIMEOUT 300)