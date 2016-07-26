
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsCopyFilesToBinaryDir)

#
# A) Define the package
#

TRIBITS_PACKAGE(Core)
SET(${PACKAGE_NAME}_VERSION "1.0")

#
# B) Set up package-specific options
#

ADD_SUBDIRECTORY(Geometry)
SET(HEADERS ${HEADERS} ${Geometry_HEADERS})
SET(SOURCES ${SOURCES} ${Geometry_SOURCES})

ADD_SUBDIRECTORY(HyperReduction)
SET(HEADERS ${HEADERS} ${HyperReduction_HEADERS})
SET(SOURCES ${SOURCES} ${HyperReduction_SOURCES})

ADD_SUBDIRECTORY(ReducedBasis)
SET(HEADERS ${HEADERS} ${ReducedBasis_HEADERS})
SET(SOURCES ${SOURCES} ${ReducedBasis_SOURCES})

ADD_SUBDIRECTORY(Utilities)
SET(HEADERS ${HEADERS} ${Utilities_HEADERS})
SET(SOURCES ${SOURCES} ${Utilities_SOURCES})

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
  tauruscore
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)

TRIBITS_ADD_TEST_DIRECTORIES(Testsuite)

#
# Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
