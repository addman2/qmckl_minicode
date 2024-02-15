add_test(NAME "QMCkl GPU reader"
  COMMAND $<TARGET_FILE:qmckl_gpu_example_trexio> h.hdf5 qmckl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
