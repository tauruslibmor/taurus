SET(HyperReduction_HEADERS
  HyperReduction/ImportRhsSnapshots.hpp
  HyperReduction/ImportMatrixSnapshots.hpp
  HyperReduction/DEIM.hpp
  HyperReduction/MeshConnectivities.hpp
  HyperReduction/VectorizedToMatrices.hpp
  HyperReduction/ReducedMesh.hpp
  HyperReduction/HyperReductionVector.hpp
  HyperReduction/HyperReductionMatrix.hpp
  HyperReduction/BuildReducedMesh.hpp
CACHE INTERNAL "")

SET(HyperReduction_SOURCES
  HyperReduction/ImportRhsSnapshots.cpp
  HyperReduction/ImportMatrixSnapshots.cpp
  HyperReduction/MeshConnectivities.cpp
  HyperReduction/VectorizedToMatrices.cpp
  HyperReduction/HyperReductionVector.cpp
  HyperReduction/HyperReductionMatrix.cpp
  HyperReduction/BuildReducedMesh.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${HyperReduction_HEADERS})
  STRING(REPLACE "HyperReduction/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/taurus/Core/HyperReduction"
  COMPONENT ${PACKAGE_NAME}
)
