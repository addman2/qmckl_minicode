add_test(NAME "QMCkl reader"
  COMMAND $<TARGET_FILE:qmckl_reader> h.hdf5 qmckl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
