cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

project(TurboRVB_minicode LANGUAGES Fortran C) 
find_package (Python COMPONENTS Interpreter)

include(CTest)

set(OFFLOAD_FLAGS   ""    CACHE STRING    "Flags of openmp offload")

find_package(OpenMP COMPONENTS Fortran)
if (${OpenMP_Fortran_FOUND})
  message(STATUS "OpenMP offload flags: ${OFFLOAD_FLAGS}")
  if (NOT "${OFFLOAD_FLAGS}" STREQUAL "")
    set(OpenMP_Fortran_FLAGS ${OFFLOAD_FLAGS})
  endif()
endif()

message(STATUS "OpenMP Fortran flags: ${OpenMP_Fortran_FLAGS}")

# Set preprocessor directives
set(PPDIRECTIVES "")

find_package(TREXIO)
if (TREXIO_FOUND)
    list( APPEND PPDIRECTIVES _TREXIO )
endif()

find_package(QMCKL)
if (QMCKL_FOUND)
    list( APPEND PPDIRECTIVES _QMCKL )
endif()

find_package(QMCKL_GPU)
if (QMCKL_GPU_FOUND)
    list( APPEND PPDIRECTIVES _QMCKL_GPU )
endif()

# Make special target that will create src/makefun.F90 file
add_custom_command(
    OUTPUT makefun.F90
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/assemble/
    COMMAND ${CMAKE_COMMAND} -E env ${Python_EXECUTABLE} assemble.py > ${CMAKE_CURRENT_BINARY_DIR}/makefun.F90
    COMMENT "Generating makefun.F90"
    VERBATIM
    DEPENDS ${CMAKE_SOURCE_DIR}/src/assemble/orb_*.f90 ${CMAKE_SOURCE_DIR}/src/assemble/assemble.py
    )

list(APPEND EXECUTABLES qmckl_reader)
list(APPEND EXECUTABLES comparator)

add_library(turborvb_minicode_lib STATIC
    src/provider_qmckl.F90
    src/provider_qmckl_gpu.F90
    src/provider_makefun.F90
    src/constants.F90
    src/qmckl.F90
    src/qmckl_gpu.F90
    src/trexio.F90
    src/qmckl_helpers.F90
    src/qmckl_gpu_helpers.F90
    ${CMAKE_CURRENT_BINARY_DIR}/makefun.F90
    )

foreach(EXECUTABLE ${EXECUTABLES})
    add_executable(${EXECUTABLE})
    target_sources(${EXECUTABLE}
        PRIVATE
        src/${EXECUTABLE}.F90
        )
    target_link_libraries(${EXECUTABLE} PRIVATE turborvb_minicode_lib)
endforeach()

message(STATUS "Preprocessor directives: ${PPDIRECTIVES}")

list(APPEND EXECUTABLES turborvb_minicode_lib)

foreach(EXECUTABLE IN LISTS EXECUTABLES)
  set_target_properties(${EXECUTABLE} PROPERTIES LINK_FLAGS "-g -O0 ${OpenMP_Fortran_FLAGS}")

  if (OpenMP_Fortran_FOUND)
    #target_link_libraries(${EXECUTABLE} PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  if (TREXIO_FOUND)
    target_link_libraries(${EXECUTABLE} PRIVATE ${TREXIO_LIBRARIES})
    target_include_directories(${EXECUTABLE} PRIVATE ${TREXIO_INCLUDE_DIRS})
  endif()

  if (QMCKL_FOUND)
    target_link_libraries(${EXECUTABLE} PRIVATE ${QMCKL_LIBRARIES})
    target_include_directories(${EXECUTABLE} PRIVATE ${QMCKL_INCLUDE_DIRS})
  endif()

  if (QMCKL_FOUND)
    target_link_libraries(${EXECUTABLE} PRIVATE ${QMCKL_GPU_LIBRARIES})
    target_include_directories(${EXECUTABLE} PRIVATE ${QMCKL_GPU_INCLUDE_DIRS})
  endif()

  get_target_property(SRC_FILES ${EXECUTABLE} SOURCES)
  foreach(SRC_FILE ${SRC_FILES})
    foreach(PPDIRECTIVE ${PPDIRECTIVES})
      set_property(SOURCE ${SRC_FILE} APPEND PROPERTY COMPILE_DEFINITIONS ${PPDIRECTIVE})
      message(STATUS "Adding preprocessor directive ${PPDIRECTIVE} to ${SRC_FILE}")
    endforeach()
    set_source_files_properties(${SRC_FILE} PROPERTIES Fortran_PREPROCESS ON)
    set_source_files_properties(${SRC_FILE} PROPERTIES Fortran_FREE ON)
    # Flags for checking
    # Set OpenMP flags for Fortran files, supporting GPU acceleration
    set (COMPILE_FLAGS "-g -O0")
    if (OpenMP_Fortran_FOUND)
      message(STATUS "Adding OpenMP flags to ${SRC_FILE}: ${OpenMP_Fortran_FLAGS}")
      set(COMPILE_FLAGS "${COMPILE_FLAGS} ${OpenMP_Fortran_FLAGS}")
    endif()
    set_source_files_properties(${SRC_FILE} PROPERTIES COMPILE_FLAGS ${COMPILE_FLAGS})
  endforeach()
endforeach()

# Add tests ##############################################################

enable_testing()

#add_subdirectory(tests/with_qmckl_h)
#add_subdirectory(tests/makefun)
add_subdirectory(tests/reader)
add_subdirectory(tests/comparator)
#add_subdirectory(tests/gpu_cpu_compare)

if (NOT ${OFFLOAD_FLAGS} STREQUAL "")
  #add_subdirectory(tests/openmp)
endif()

