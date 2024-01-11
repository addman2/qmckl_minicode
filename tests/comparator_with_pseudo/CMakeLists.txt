add_test(NAME "ComparatorWithPseudo"
  COMMAND $<TARGET_FILE:comparator_with_pseudo>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
