# File for module specific CMake tests.
## my code for o2scl  ######################
find_package(o2scl REQUIRED)
include_directories(${o2scl_INCLUDE_DIRS})
set(LIBS ${LIBS} ${o2scl_LIBRARIES})
############################################

