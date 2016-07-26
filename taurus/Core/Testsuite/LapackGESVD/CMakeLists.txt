INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE(
  Test_GESVD
  SOURCES main.cpp
  COMM serial mpi
)
