
find_path(o2scl_INCLUDE_DIR o2scl/exception.h)
find_library(o2scl_LIBRARY o2scl)

set(o2scl_LIBRARIES ${o2scl_LIBRARY} )
set(o2scl_INCLUDE_DIRS ${o2scl_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(o2scl  DEFAULT_MSG
                                  o2scl_LIBRARY o2scl_INCLUDE_DIR)

				
				message(">>>>>>>>>>>> ${o2scl_LIBRARIES}")
				message(">>>>>>>>>>>> ${o2scl_INCLUDE_DIR}")
				
mark_as_advanced(o2scl_INCLUDE_DIR o2scl_LIBRARY )
