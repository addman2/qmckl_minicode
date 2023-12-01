add_test(NAME "With QMCkl (hydrogen)"
  COMMAND $<TARGET_FILE:turborvb_minicode> h.hdf5 qmckl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
