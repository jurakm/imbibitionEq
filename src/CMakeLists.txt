
add_executable("imbibitionEq" imbibitionEq.cc)

add_library("libImbibition"  SHARED  imbibitionFunctions.hh integration.hh linear_analytic.hh mgf.hh parameters.hh tableXY.hh 
             text_output.hh timeMng.hh vangenuchten.hh)

target_link_libraries("libImbibition"  "${o2scl_LIBRARIES}")
SET_TARGET_PROPERTIES("libImbibition" PROPERTIES LINKER_LANGUAGE CXX)

set_property(TARGET imbibitionEq APPEND_STRING PROPERTY COMPILE_FLAGS "-I${o2scl_INCLUDE_DIRS}")
#set_property(TARGET imbibitionEq APPEND PROPERTY LINK_LIBRARIES "${o2scl_LIBRARIES}")
set_property(TARGET imbibitionEq APPEND PROPERTY LINK_LIBRARIES "libImbibition")

####  my code ############################
#target_link_libraries("imbibitionEq" ${LIBS})
##########################################
target_link_dune_default_libraries("imbibitionEq")
