SET(ReducedBasis_HEADERS
  ReducedBasis/ImportSnapshots.hpp
  ReducedBasis/V_POD.hpp
  ReducedBasis/BuildReducedBasis.hpp
  ReducedBasis/ReductionVector.hpp
  ReducedBasis/BlockDiagonalBasis.hpp
CACHE INTERNAL "")

SET(ReducedBasis_SOURCES
  ReducedBasis/ImportSnapshots.cpp
  ReducedBasis/BuildReducedBasis.cpp
  ReducedBasis/ReductionVector.cpp
  ReducedBasis/BlockDiagonalBasis.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${ReducedBasis_HEADERS})
  STRING(REPLACE "ReducedBasis/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/taurus/Core/ReducedBasis"
  COMPONENT ${PACKAGE_NAME}
)
