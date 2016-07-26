SET(Models_HEADERS
  Models/StructuralModel.hpp
  Models/LinearElasticMaterial.hpp
  Models/StVenantKirchhoffMaterial.hpp
  Models/AssemblyElementalStructure.hpp
  Models/ExpressionDefinitions.hpp
  Models/StructuralSolverTierI.hpp
  Models/StructuralSolverTierII.hpp
  Models/StructuralSolverTierIII.hpp
CACHE INTERNAL "")

SET(Models_SOURCES
  Models/StructuralModel.cpp
  Models/LinearElasticMaterial.cpp
  Models/StVenantKirchhoffMaterial.cpp
  Models/AssemblyElementalStructure.cpp
  Models/ExpressionDefinitions.cpp
  Models/StructuralSolverTierI.cpp
  Models/StructuralSolverTierII.cpp
  Models/StructuralSolverTierIII.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Models_HEADERS})
  STRING(REPLACE "Models/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

#INSTALL(
#  FILES ${LOCAL_HEADERS}
#  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/taurus/CSM/Models"
#  COMPONENT ${PACKAGE_NAME}
#)
