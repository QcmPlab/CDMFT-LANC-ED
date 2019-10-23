CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(dmft_ed Fortran)
SET(VERSION 0.15)


SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Check that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
OPTION(USE_MPI "Set MPI pre-compilation flag" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
OPTION(WITH_BLAS_LAPACK "Use internal Blas/Lapack" OFF)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)



############################################################
# Look for MPI if needed
############################################################
# reSet compiler to MPI is required
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF(MPI_Fortran_FOUND)
    SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    MESSAGE(STATUS "Set Fortran compiler to ${Yellow}${CMAKE_Fortran_COMPILER}${ColourReset}, ID=${CMAKE_Fortran_COMPILER_ID}")
    SET(MPI_CPP "MPI")		#pre-processing option
  ELSE()
    MESSAGE(FATAL_ERROR "MPI Found but No MPI-Fortran compiler can be determined.")    
  ENDIF()
ELSE(USE_MPI)
  SET(MPI_CPP "")
ENDIF(USE_MPI)




############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_SHA1
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "Git version: ${GIT_SHA1}" )
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_BRANCH
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    SET(GIT_BRANCH_ "_${GIT_BRANCH}")
    MESSAGE( STATUS "Git branch: ${GIT_BRANCH}" )
    IF( "${GIT_BRANCH}" MATCHES "master")
      UNSET(GIT_BRANCH)
      UNSET(GIT_BRANCH_)
    ENDIF()
  ELSE(GIT_FOUND)
    SET(GIT_SHA1 0)
  ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
SET(REVISION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/revision.inc)
FILE(WRITE  ${REVISION_FILE} "character(len=41),parameter,public :: revision = \"${GIT_SHA1}\"\n")




############################################################
# Look for required libraries: blas,lapack,scifor,dmft_tools
############################################################
# Find Scifor
FIND_PACKAGE(PkgConfig)
PKG_SEARCH_MODULE(SCIFOR REQUIRED scifor)
IF(${SCIFOR_FOUND})
  MESSAGE(STATUS "SCIFOR found at: ${SCIFOR_PREFIX}")
ELSE()
  MESSAGE(FATAL_ERROR "
Can not find SciFortran in this system. Check SciFor installation. Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_scifor_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
IF("${SCIFOR_LDFLAGS}" MATCHES "mkl")
  FIND_PACKAGE(MKL QUIET)
ELSEIF("${SCIFOR_LDFLAGS}" MATCHES "blas" OR "lapack")
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})


# Find Dmft_tools
PKG_SEARCH_MODULE(DMFT_TOOLS REQUIRED dmft_tools)
IF(${DMFT_TOOLS_FOUND})
  MESSAGE(STATUS "DMFT_TOOLS found at: ${DMFT_TOOLS_PREFIX}")
ELSE()
  MESSAGE( FATAL_ERROR "
Can not find DMFT_TOOLS in this system. Check DMFT_TOOLS installation. Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_dmft_tools_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${DMFT_TOOLS_INCLUDE_DIRS})





############################################################
# Set Fortran options based on BUILD_TYPE and FC ID
############################################################
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
  SET(CMAKE_Fortran_MODDIR_FLAG   "-J")
  SET(CMAKE_Fortran_FLAGS         "-cpp -ffree-line-length-none ") # -mcmodel=large  this is to remove the 2Gb limit of virtual memory allocation
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -funroll-loops")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -fbounds-check")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3  -funroll-loops")  ##-fno-protect-parens -flto -ffast-math -march=native
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ") #remember the ending white space here 
  SET(CMAKE_Fortran_FLAGS         "-fpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -ftz")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-p -O0 -g -fpe1 -warn -debug extended -traceback -check all,noarg_temp_created")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftz")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES G95)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-fmod=")
  SET(CMAKE_Fortran_FLAGS         "-cpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O1  -fsloppy-char")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ftrace=full")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ")
  SET(CMAKE_Fortran_FLAGS         "-")
  SET(CMAKE_Fortran_FLAGS         "")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mcmodel=medium -fast -Munroll")
endif()


#set MPI preprocessing defition:
IF(USE_MPI)
  ADD_DEFINITIONS(-D_MPI)
ELSE(USE_MPI)
  ADD_DEFINITIONS(-D_)
ENDIF(USE_MPI)

IF( "${BUILD_TYPE}" MATCHES "DEBUG")
  ADD_DEFINITIONS(-D_DEBUG)
ENDIF()




############################################################
# Define the actual files and folders that make up the build
############################################################
SET(USER_HOME $ENV{HOME})

# Define some directories
SET(DRIVERS ${CMAKE_SOURCE_DIR}/drivers)

#define PREFIX
SET(PREFIX  "$ENV{HOME}/.bin" CACHE PATH "Prefix prepended to install directories")
SET(CMAKE_INSTALL_PREFIX "${PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

# define the executable to compiler
SET(EXE edn_hm_bethe CACHE STRING "Driver code to compile")
# SET(CMAKE_INSTALL_EXE "${EXE}" CACHE INTERNAL "Driver code to compile" FORCE)


FILE(GLOB DMFT_ED_LOCAL_SRC
  ED_SPARSE_MATRIX.f90
  ED_INPUT_VARS.f90
  ED_VARS_GLOBAL.f90
  ED_AUX_FUNX.f90
  ED_HLOC_DECOMPOSITION.f90
  ED_SETUP.f90
  ED_EIGENSPACE.f90
  ED_BATH.f90
  ED_BATH_FUNCTIONS.f90
  ED_IO.f90
  ED_HAMILTONIAN_COMMON.f90
  ED_HAMILTONIAN_SPARSE_HxV.f90
  ED_HAMILTONIAN_DIRECT_HxV.f90
  ED_HAMILTONIAN.f90
  ED_DIAG.f90
  ED_GF_SHARED.f90
  ED_GF_NORMAL.f90
  #ED_GF_CHISPIN.f90
  ED_GREENS_FUNCTIONS.f90
  ED_OBSERVABLES.f90    
  ED_FIT_CHI2.f90
  ED_MAIN.f90
  CDMFT_ED.f90
  )

SET(DMFT_ED_OBJS ${DMFT_ED_LOCAL_SRC})

ADD_EXECUTABLE(dmft_ed ${DRIVERS}/${EXE}.f90 ${DMFT_ED_OBJS})


target_LINK_LIBRARIES(dmft_ed ${DMFT_TOOLS_LDFLAGS} ${SCIFOR_LDFLAGS})
SET_TARGET_PROPERTIES(dmft_ed PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}" OUTPUT_NAME "${EXE}${GIT_BRANCH_}")
MESSAGE(STATUS "Executable created at: ${Yellow}${CMAKE_INSTALL_PREFIX}/${EXE}${GIT_BRANCH_}${ColourReset}")

ADD_CUSTOM_TARGET(
  completion
  COMMAND scifor_completion.sh ${DRIVERS}/${EXE}.f90
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
