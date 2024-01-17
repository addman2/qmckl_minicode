add_test(NAME "MakefunTester"
  COMMAND $<TARGET_FILE:makefun_tester>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
