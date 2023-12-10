add_test(NAME "With QMCkl (hydrogen)"
  COMMAND $<TARGET_FILE:makefun_runner> h.hdf5 qmckl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
