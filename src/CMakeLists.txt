
add_executable("imbibitionEq" imbibitionEq.cc)

set_property(TARGET imbibitionEq APPEND_STRING PROPERTY COMPILE_FLAGS "-I${o2scl_INCLUDE_DIRS}")
set_property(TARGET imbibitionEq APPEND PROPERTY LINK_LIBRARIES "${o2scl_LIBRARIES}")

####  my code ############################
#target_link_libraries("imbibitionEq" ${LIBS})
##########################################
target_link_dune_default_libraries("imbibitionEq")
