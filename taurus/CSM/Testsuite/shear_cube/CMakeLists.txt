INCLUDE(TribitsAddExecutable)
INCLUDE(AddSubdirectories)

TRIBITS_ADD_EXECUTABLE(
  GenerateSolutionSnapshots_shear_cube
  SOURCES main_GenerateSolutionSnapshots.cpp
  COMM serial mpi
) 

TRIBITS_ADD_EXECUTABLE(
  GenerateSystemSnapshots_shear_cube
  SOURCES main_GenerateSystemSnapshots.cpp
  COMM serial mpi
) 


TRIBITS_ADD_EXECUTABLE(
  GenerateSolutionBasis_shear_cube
  SOURCES main_GenerateSolutionBasis.cpp
  COMM serial mpi
) 

TRIBITS_ADD_EXECUTABLE(
  GenerateHyperROM_shear_cube
  SOURCES main_GenerateHyperROM.cpp
  COMM serial mpi
) 

TRIBITS_ADD_EXECUTABLE(
OnlineCompare_shear_cube
SOURCES main_OnlineCompare.cpp
COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
SolveOnline_shear_cube
SOURCES main_SolveOnline.cpp
COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  launcher_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES launch_all.sh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  LHS_ParametersConfigurations
  CREATE_SYMLINK
  SOURCE_FILES LHS_ParametersConfigurations.m
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Structure_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Structure_shear_cube_mesh
  CREATE_SYMLINK
  SOURCE_FILES Cube.mesh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Structure_shear_cube_msh
  CREATE_SYMLINK
  SOURCE_FILES Cube.msh
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Structure_shear_cube_geo
  CREATE_SYMLINK
  SOURCE_FILES Cube.geo
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Hyper_Structure_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES SolverParamList.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Structure_online_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES data_online
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Structure_onlineCompare_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES data_compare
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Structure_TP_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES TrainingParameters.vit
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  VisualizationParameters_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES VisualizationParameters.vit
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  TestingParameters_shear_cube
  CREATE_SYMLINK
  SOURCE_FILES TestingParameters.vit
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
