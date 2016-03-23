#!/bin/bash

echo "Running SDSS_CF_Parallel..." 
time ./SDSS_CF_Parallel

echo "Copying Files to DMCorrelationFunctionParallel"
cp output_Xi*.txt ../DMCorrelationFunctionParallel

echo "Running Plotting Software for SDSS plots..."
time python Plot_Correlation_Func_From_C.py

echo "Running DM_CF_Parallel..."
time ../DMCorrelationFunctionParallel/DM_CF_Parallel

echo "Running Plotting Software for DM plots..."
time python ../DMCorrelationFunctionParallel/Plot_Correlation_Func_from_C2.py

echo "Copying plots to HOME directory."
cp ../DMCorrelationFunctionParallel/*.png ~

echo "FINISHED RUN_FULL_CF.sh" 






