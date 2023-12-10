add_test(NAME "With QMCkl cpu - gpu (hydrogen)"
  COMMAND $<TARGET_FILE:qmckl_gpu_cpu_compare> h.hdf5
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
