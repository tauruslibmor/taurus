#!/bin/bash

# usage: sh launch_all num_processor

mpirun -n $1 ./CSM_GenerateSolutionSnapshots_shear_cube.exe -f data | tee output_GenerateSolutionSnapshots.txt
mpirun -n $1 ./CSM_GenerateSolutionBasis_shear_cube.exe -f data | tee output_GenerateSolutionBasis.txt
mpirun -n $1 ./CSM_GenerateSystemSnapshots_shear_cube.exe -f data | tee output_GenerateSystemSnapshots.txt
mpirun -n $1 ./CSM_GenerateHyperROM_shear_cube.exe -f data | tee output_GenerateHyperROM.txt
mpirun -n $1 ./CSM_SolveOnline_shear_cube.exe -f data_online | tee output_SolveOnline.txt
mpirun -n $1 ./CSM_OnlineCompare_shear_cube.exe -f data_compare | tee output_OnlineCompare.txt
