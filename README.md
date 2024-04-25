Here's a pipeline to treat downloaded SRA files from NCBI to get in the end a count of each allele per position.

The scripts are made to run on HPC servers using the SLURM environment. 
To make this work, you can either launch the Treat_single_library.sh bash script alone or specify a list of libraries to the Treat_List_Of_library.sh bash script. 
The complete path of the SRA file need to be present in the text file that will be given, including the ".sra suffix". The number of parallel jobs can be changed manually in the Treat_List_Of_library.sh file.


Here are the current dependencies :

SRA toolkit (v 2.1.0.8)
BWA v0.7.17
samtools v1.17
htslib v1.17

TrimGalore! v0.6.0
CutAdapt v4.2

iVar v1.3
pileup2base perl script (https://github.com/riverlee/pileup2base)

Qualimap v2.2.1

The pipeline could also work with other program versions.
