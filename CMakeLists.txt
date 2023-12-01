cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

project(TurboRVB_minicode LANGUAGES Fortran C) 

include(CTest)

set(OFFLOAD_FLAGS   ""    CACHE STRING    "Flags of openmp offload")

add_executable(turborvb_minicode)

add_library(turborvb_minicode_lib STATIC
    src/constants.F90
    src/qmckl.F90
    src/trexio.F90
    )

target_sources(turborvb_minicode
    PRIVATE
    src/makefun.F90
    src/test_makefun_qmckl.F90
    )

# Set preprocessor directives
set(PPDIRECTIVES "")

find_package(OpenMP COMPONENTS Fortran)
if (${OpenMP_Fortran_FOUND})
  if (NOT "${OFFLOAD_FLAGS}" STREQUAL "")
    target_link_libraries(turborvb_minicode_lib PUBLIC ${OFFLOAD_FLAGS})
    target_link_libraries(turborvb_minicode PUBLIC ${OFFLOAD_FLAGS})
    set_property(TARGET turborvb_minicode_lib APPEND PROPERTY LINK_FLAGS ${OFFLOAD_FLAGS})
    set_property(TARGET turborvb_minicode APPEND PROPERTY LINK_FLAGS ${OFFLOAD_FLAGS})
    set(OpenMP_Fortran_FLAGS ${OFFLOAD_FLAGS})
  else()
    target_link_libraries(turborvb_minicode_lib PUBLIC OpenMP::OpenMP_Fortran)
    target_link_libraries(turborvb_minicodeb PUBLIC OpenMP::OpenMP_Fortran)
  endif()
endif()

find_package(TREXIO)
if (TREXIO_FOUND)
    list( APPEND PPDIRECTIVES _TREXIO )
    target_link_libraries(turborvb_minicode_lib PUBLIC ${TREXIO_LIBRARIES})
    include_directories( ${TREXIO_INCLUDE_DIRS} )
endif()

find_package(QMCKL)
if (QMCKL_FOUND)
    list( APPEND PPDIRECTIVES _QMCKL )
    target_link_libraries(turborvb_minicode_lib PUBLIC ${QMCKL_LIBRARIES})
    include_directories( ${QMCKL_INCLUDE_DIRS} )
endif()

#find_package(QMCKL_GPU)
#if (QMCKL_GPU_FOUND)
#    list( APPEND PPDIRECTIVES _QMCKL_GPU )
#    link_libraries(${QMCKL_GPU_LIBRARIES})
#    include_directories( ${QMCKL_GPU_INCLUDE_DIRS} )
#endif()

message(STATUS "Preprocessor directives: ${PPDIRECTIVES}")

get_target_property(SRC_FILES turborvb_minicode_lib SOURCES)
foreach(SRC_FILE ${SRC_FILES})
  foreach(PPDIRECTIVE ${PPDIRECTIVES})
    set_property(SOURCE ${SRC_FILE} APPEND PROPERTY COMPILE_DEFINITIONS ${PPDIRECTIVE})
  endforeach()
  set_source_files_properties(${SRC_FILE} PROPERTIES Fortran_PREPROCESS ON)
  set_source_files_properties(${SRC_FILE} PROPERTIES Fortran_FREE ON)
  # Set OpenMP flags for Fortran files, supporting GPU acceleration
  if (OpenMP_Fortran_FOUND)
    set_source_files_properties(${SRC_FILE} PROPERTIES COMPILE_FLAGS ${OpenMP_Fortran_FLAGS})
  endif()
endforeach()

get_target_property(SRC_FILES turborvb_minicode SOURCES)
foreach(SRC_FILE ${SRC_FILES})
  foreach(PPDIRECTIVE ${PPDIRECTIVES})
    set_property(SOURCE ${SRC_FILE} APPEND PROPERTY COMPILE_DEFINITIONS ${PPDIRECTIVE})
  endforeach()
  set_source_files_properties(${SRC_FILE} PROPERTIES Fortran_PREPROCESS ON)
  set_source_files_properties(${SRC_FILE} PROPERTIES Fortran_FREE ON)
  if (OpenMP_Fortran_FOUND)
    set_source_files_properties(${SRC_FILE} PROPERTIES COMPILE_FLAGS ${OpenMP_Fortran_FLAGS})
  endif()
endforeach()

# Link libraries
target_link_libraries(turborvb_minicode PRIVATE turborvb_minicode_lib)

add_subdirectory(tests/with_qmckl_h)
if (NOT ${OFFLOAD_FLAGS} STREQUAL "")
  add_subdirectory(tests/openmp)
endif()
