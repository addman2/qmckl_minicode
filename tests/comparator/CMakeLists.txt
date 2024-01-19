add_test(NAME "ComparatorWithoutPseudo"
  COMMAND $<TARGET_FILE:comparator>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
