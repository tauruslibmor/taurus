SET(Utilities_HEADERS
  Utilities/DOF_Extractor.hpp
  Utilities/Utilities.hpp
  Utilities/DenseHDF5.hpp
  Utilities/VITreader.hpp
  Utilities/assembly_functor.hpp
  Utilities/GPmodel.hpp
  Utilities/Timer.hpp
  Utilities/TimeAndExtrapolationHandler.hpp
CACHE INTERNAL "")

SET(Utilities_SOURCES
  Utilities/DOF_Extractor.cpp
  Utilities/Utilities.cpp
  Utilities/DenseHDF5.cpp
  Utilities/VITreader.cpp
  Utilities/GPmodel.cpp
  Utilities/Timer.cpp
  Utilities/TimeAndExtrapolationHandler.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Utilities_HEADERS})
  STRING(REPLACE "Utilities/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/taurus/Core/Utilities"
  COMPONENT ${PACKAGE_NAME}
)
