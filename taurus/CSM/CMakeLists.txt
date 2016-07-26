
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

TRIBITS_PACKAGE(CSM)
SET(${PACKAGE_NAME}_VERSION "1.0")

#
# B) Set up package-specific options
#

ADD_SUBDIRECTORY(Models)
SET(HEADERS ${HEADERS} ${Models_HEADERS})
SET(SOURCES ${SOURCES} ${Models_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  tauruscsm
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)

TRIBITS_ADD_TEST_DIRECTORIES(Testsuite)

#
# Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
