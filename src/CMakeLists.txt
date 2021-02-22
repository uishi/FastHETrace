cmake_minimum_required (VERSION 3.5.1)
                                                                                                       
project(demo CXX)
set(CMAKE_CXX_STANDARD 11)

option(USE_FAST_AUTO "Set to ON to use table-based automorphism in PALISADE" ON)

#find_package(Palisade 1.9.2 EXACT REQUIRED PATHS $ENV{HOME}/lib/CMake/)
find_package(Palisade 1.9.2 EXACT REQUIRED)
if (${Palisade_FOUND})
  message("palisade found")
  message( ${PALISADE_INCLUDE})
  message( ${PALISADE_LIBDIR})
else (${Palisade_FOUND})
  message(FATAL_ERROR "palisade not found")
endif (${Palisade_FOUND})

message( ${PALISADE_INCLUDE})
message( ${PALISADE_INCLUDE}/core )
message( ${PALISADE_INCLUDE}/pke )


set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )

include_directories(".")
include_directories("../include")
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )

link_directories( ${PALISADE_LIBDIR} )
link_libraries( ${PALISADE_LIBRARIES} )


set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )

if (${USE_FAST_AUTO})
  message(" Fast automorph option ON ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FAST_AUTOMORPH")
endif()

include_directories( ${OPENMP_INCLUDES} )
link_directories( ${OPENMP_LIBRARIES} )
message("CMAKE FLAGS = ${CMAKE_CXX_FLAGS}")

add_executable(benchunroll_dynamic_dnum  benchunroll_dynamic_dnum.cpp)
add_executable(bench_seq_dynamic_dnum  bench_seq_dynamic_dnum.cpp)
