add_test(NAME "OpenMP data offload test"
         COMMAND $<TARGET_FILE:turborvb_minicode> h.hdf5 openmp
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
         )
