INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  GP_model_Test
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  training_set_GP_model_Test
  CREATE_SYMLINK
  SOURCE_FILES GPtraining.txt
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
