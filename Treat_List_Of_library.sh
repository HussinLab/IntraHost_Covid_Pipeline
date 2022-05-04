#!/bin/bash

if [ $CC_CLUSTER == "graham" ]; then 
	parallel --jobs 32 /project/6005588/shared/bin/covid19-dev/IntraHost_Covid_Pipeline/Treat_single_library.sh {1} :::: $1 
elif [ $CC_CLUSTER == "narval" ]; then 
	parallel --jobs 32 /lustre06/project/6065672/shared/covid-19/RNAseq/viral/NCBI/IntraHost_Covid_Pipeline/Treat_single_library.sh {1} :::: $1 
fi
