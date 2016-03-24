#!/bin/bash

echo "Deleting previous output files!"
rm ./output*.txt 
rm ../DMCorrelationFunctionParallel/output*.txt

echo "Setting the Number of Threads for PRAGMA *OOOMP*."
export OMP_NUM_THREADS=32 

echo "Running SDSS_CF_Parallel..." 
time ./SDSS_CF_Parallel

echo "Copying Files to DMCorrelationFunctionParallel for Bias Function"
cp output_Xi*.txt ../DMCorrelationFunctionParallel

echo "Running DM_CF_Parallel..."
time ../DMCorrelationFunctionParallel/DM_CF_Parallel

echo "FINISHED CALCULATIONS!"

echo "Running Plotting Software for SDSS plots..."
time python Plot_Correlation_Func_From_C.py

echo "Running Plotting Software for DM plots..."
time python ../DMCorrelationFunctionParallel/Plot_Correlation_Func_From_C2.py

echo "Copying SDSS and DM plots to HOME directory. "
cp ../DMCorrelationFunctionParallel/HW3*Para*.png ~
cp ./HW3*Para*png ~

echo "FINISHED RUN_FULL_CF.sh" 






