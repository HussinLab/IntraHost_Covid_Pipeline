#!/bin/bash
parallel --jobs 32 /lustre06/project/6065672/shared/covid-19/RNAseq/viral/NCBI/IntraHost_Covid_Pipeline/Treat_single_library.sh {1} :::: $1
