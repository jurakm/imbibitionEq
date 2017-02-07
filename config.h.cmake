/* config.h.  Generated from config.h.cmake by CMake.
   It was generated from config.h.cmake which in turn is generated automatically
   from the config.h.cmake files of modules this module depends on. */

/* Define to 1 if you have module imbibitionEq available */
#cmakedefine01 HAVE_IMBIBITIONEQ


/* Define to 1 if you have module dune-common available */
#cmakedefine01 HAVE_DUNE_COMMON


/* Define to 1 if you have module dune-geometry available */
#cmakedefine01 HAVE_DUNE_GEOMETRY


/* Define to 1 if you have module dune-typetree available */
#cmakedefine01 HAVE_DUNE_TYPETREE


/* Define to 1 if you have module dune-istl available */
#cmakedefine01 HAVE_DUNE_ISTL


/* Define to 1 if you have module dune-localfunctions available */
#cmakedefine01 HAVE_DUNE_LOCALFUNCTIONS


/* Define to 1 if you have module dune-grid available */
#cmakedefine01 HAVE_DUNE_GRID


/* Define to 1 if you have module dune-pdelab available */
#cmakedefine01 HAVE_DUNE_PDELAB


/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "${DUNE_COMMON_VERSION}"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR ${DUNE_COMMON_VERSION_MAJOR}

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR ${DUNE_COMMON_VERSION_MINOR}

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION ${DUNE_COMMON_VERSION_REVISION}

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL ${DUNE_MINIMAL_DEBUG_LEVEL}

/* does the compiler support __attribute__((deprecated))? */
#cmakedefine HAS_ATTRIBUTE_DEPRECATED 1

/* does the compiler support __attribute__((deprecated("message"))? */
#cmakedefine HAS_ATTRIBUTE_DEPRECATED_MSG 1

/* does the compiler support __attribute__((unused))? */
#cmakedefine HAS_ATTRIBUTE_UNUSED 1

/* old feature support macros which were tested until 2.3, keep around for one more release */
/* As these are now always supported due to the new compiler requirements, they are directly */
/* defined without an explicit test. */
#define HAVE_INTEGRAL_CONSTANT 1
#define HAVE_DUNE_HASH 1
#define HAVE_STD_HASH 1
#define HAVE_TYPE_TRAITS 1
#define HAVE_VARIADIC_TEMPLATES  1
#define HAVE_VARIADIC_CONSTRUCTOR_SFINAE 1
#define HAVE_RVALUE_REFERENCES 1

/* Define if you have a BLAS library. */
#cmakedefine HAVE_BLAS 1

/* Define to ENABLE_BOOST if the Boost library is available */
#cmakedefine HAVE_DUNE_BOOST ENABLE_BOOST

/* does the compiler support abi::__cxa_demangle */
#cmakedefine HAVE_CXA_DEMANGLE 1

/* Define if you have LAPACK library. */
#cmakedefine HAVE_LAPACK 1

/* Define to 1 if you have the <malloc.h> header file. */
// Not used! #cmakedefine01 HAVE_MALLOC_H

/* Define if you have the MPI library. This is only true if MPI was found by
   configure _and_ if the application uses the DUNEMPICPPFLAGS (or the
   deprecated MPI_CPPFLAGS) */
#cmakedefine HAVE_MPI ENABLE_MPI

/* Define if you have the GNU GMP library. The value should be ENABLE_GMP
   to facilitate activating and deactivating GMP using compile flags. */
#cmakedefine HAVE_GMP ENABLE_GMP

/* Define to 1 if you have the symbol mprotect. */
#cmakedefine HAVE_MPROTECT 1

/* Define to 1 if nullptr is supported */
#cmakedefine HAVE_NULLPTR 1

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H 1

/* Define to 1 if you have <sys/mman.h>. */
#cmakedefine HAVE_SYS_MMAN_H 1

/* Define to 1 if the MPI2 Standard is supported */
#cmakedefine MPI_2 1



/* Define to 1 if C++11 constexpr is supported */
#cmakedefine HAVE_CONSTEXPR 1

/* Define to 1 if C++11 range-based for is supported */
#cmakedefine HAVE_RANGE_BASED_FOR 1

/* Define to 1 if C++11 nonexcept specifier is supported */
#cmakedefine HAVE_NOEXCEPT_SPECIFIER 1

/* Define to 1 if C++11 std::declval() is supported */
#cmakedefine HAVE_STD_DECLVAL 1

/* Define to 1 if the compiler properly supports testing for operator[] */
#cmakedefine HAVE_IS_INDEXABLE_SUPPORT 1

/* does the compiler support the keyword 'final'? */
#cmakedefine HAVE_KEYWORD_FINAL 1

/* Define to if the UMFPack library is available */
#cmakedefine HAVE_UMFPACK ENABLE_UMFPACK

/* Define to ENABLE_PARMETIS if you have the Parmetis library.
   This is only true if MPI was found
   by configure _and_ if the application uses the PARMETIS_CPPFLAGS */
#cmakedefine HAVE_PARMETIS ENABLE_PARMETIS

/* Define to 1 if PT-Scotch is available */
#cmakedefine HAVE_PTSCOTCH 1

/* Include always useful headers */
#include "FC.h"
#define FC_FUNC FC_GLOBAL_





/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "${DUNE_GEOMETRY_VERSION}"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR ${DUNE_GEOMETRY_VERSION_MAJOR}

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR ${DUNE_GEOMETRY_VERSION_MINOR}

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION ${DUNE_GEOMETRY_VERSION_REVISION}





/* Define to 1 if std::initializer_list is supported. */
#cmakedefine HAVE_INITIALIZER_LIST 1

/* Define to 1 if template aliases are supported. */
#cmakedefine HAVE_TEMPLATE_ALIASES 1

/* Define to 1 if decltype if supported. */
#cmakedefine HAVE_STD_DECLTYPE 1

/* Define to 1 if GCC's __typeof__ extension is supported. */
#cmakedefine HAVE_GCC___TYPEOF__ 1






/* define if the Boost::Fusion headers are available */
#cmakedefine HAVE_BOOST_FUSION

/* Define to ENABLE_BOOST if the Boost library is available */
#cmakedefine HAVE_BOOST ENABLE_BOOST

/* Define to ENABLE_SUPERLU if the SuperLU library is available */
#cmakedefine HAVE_SUPERLU ENABLE_SUPERLU

/* define to 1 because older versions of SuperLU are no longer supported*/
#define SUPERLU_POST_2005_VERSION 1

/* Define to 1 if 'expansions' is a member of 'mem_usage_t'. */
#cmakedefine HAVE_MEM_USAGE_T_EXPANSIONS @HAVE_MEM_USAGE_T_EXPANSIONS@

/* define to 1 if SuperLU header slu_ddefs.h contains SLU_DOUBLE */
#cmakedefine SUPERLU_MIN_VERSION_4_3 @SUPERLU_MIN_VERSION_4_3@

/* Define to the version of dune-istl */
#define DUNE_ISTL_VERSION "${DUNE_ISTL_VERSION}"

/* Define to the major version of dune-istl */
#define DUNE_ISTL_VERSION_MAJOR ${DUNE_ISTL_VERSION_MAJOR}

/* Define to the minor version of dune-istl */
#define DUNE_ISTL_VERSION_MINOR ${DUNE_ISTL_VERSION_MINOR}

/* Define to the revision of dune-istl */
#define DUNE_ISTL_VERSION_REVISION ${DUNE_ISTL_VERSION_REVISION}





/* Define to the version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION "${DUNE_LOCALFUNCTIONS_VERSION}"

/* Define to the major version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MAJOR ${DUNE_LOCALFUNCTIONS_VERSION_MAJOR}

/* Define to the minor version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MINOR ${DUNE_LOCALFUNCTIONS_VERSION_MINOR}

/* Define to the revision of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_REVISION ${DUNE_LOCALFUNCTIONS_VERSION_REVISION}





/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "${DUNE_GRID_VERSION}"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR ${DUNE_GRID_VERSION_MAJOR}

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR ${DUNE_GRID_VERSION_MINOR}

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION ${DUNE_GRID_VERSION_REVISION}

/* If this is set, public access to the implementation of facades like Entity,
   Geometry, etc. is granted. */
#cmakedefine DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS 1

/* This is only true if grape was found by configure _and_ if the
   application uses the flags set by add_dune_grape_flags */
#cmakedefine HAVE_GRAPE ENABLE_GRAPE

/* Define to 1 if psurface library is found */
#cmakedefine HAVE_PSURFACE 1

/* Define to 1 if AmiraMesh library is found */
#cmakedefine HAVE_AMIRAMESH 1

/* The namespace prefix of the psurface library (deprecated) */
#define PSURFACE_NAMESPACE psurface::

/* Define to 1 if you have at least psurface version 2.0 */
#cmakedefine HAVE_PSURFACE_2_0 1

/* This is only true if alugrid-library was found by configure _and_ if the
   application uses the ALUGRID_CPPFLAGS */
#cmakedefine HAVE_ALUGRID ENABLE_ALUGRID

/* Define to 1 if you have the <alugrid_parallel.h> header file. */
#cmakedefine HAVE_ALUGRID_PARALLEL_H @HAVE_ALUGRID_PARALLEL_H@

/* Define to 1 if you have the <alugrid_serial.h> header file. */
#cmakedefine HAVE_ALUGRID_SERIAL_H @HAVE_ALUGRID_SERIAL_H@

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0 */
#cmakedefine DUNE_ALBERTA_VERSION @DUNE_ALBERTA_VERSION@

/* This is only true if alberta-library was found by configure _and_ if the
   application uses the ALBERTA_CPPFLAGS */
#cmakedefine HAVE_ALBERTA ENABLE_ALBERTA

/* This is only true if UG was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
#cmakedefine HAVE_UG ENABLE_UG

/* Define to 1 if you have mkstemp function */
#cmakedefine01 HAVE_MKSTEMP

/* Grid type magic for DGF parser */
@GRID_CONFIG_H_BOTTOM@




/* This is only true if PETSc was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
#cmakedefine HAVE_PETSC ENABLE_PETSC

/* Define to 1 if you have the <tr1/unordered_set> header file. */
#cmakedefine HAVE_TR1_UNORDERED_SET 1

/* Define to 1 if you have the <tr1/unordered_map> header file. */
#cmakedefine HAVE_TR1_UNORDERED_MAP 1

/* Define to 1 if you have the <unordered_map> header file. */
#cmakedefine HAVE_UNORDERED_MAP 1

/* Define to 1 if you have the <unordered_set> header file. */
#cmakedefine HAVE_UNORDERED_SET 1

/* Define to 1 if decltype if supported. */
#cmakedefine HAVE_STD_DECLTYPE 1

/* Define to 1 if GCC's __typeof__ extension is supported. */
#cmakedefine HAVE_GCC___TYPEOF__ 1



/* begin imbibitionEq
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/

/* begin private */
/* Name of package */
#define PACKAGE "@DUNE_MOD_NAME"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "@DUNE_MAINTAINER@"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@DUNE_MOD_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@DUNE_MOD_NAME@ @DUNE_MOD_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@DUNE_MOD_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@DUNE_MOD_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@DUNE_MOD_VERSION@"

/* end private */

/* Define to the version of imbibitionEq */
#define IMBIBITIONEQ_VERSION "@IMBIBITIONEQ_VERSION@"

/* Define to the major version of imbibitionEq */
#define IMBIBITIONEQ_VERSION_MAJOR @IMBIBITIONEQ_VERSION_MAJOR@

/* Define to the minor version of imbibitionEq */
#define IMBIBITIONEQ_VERSION_MINOR @IMBIBITIONEQ_VERSION_MINOR@

/* Define to the revision of imbibitionEq */
#define IMBIBITIONEQ_VERSION_REVISION @IMBIBITIONEQ_VERSION_REVISION@

/* end imbibitionEq
   Everything below here will be overwritten
*/