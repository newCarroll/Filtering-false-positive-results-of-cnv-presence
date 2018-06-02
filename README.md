# Filtering-false-positive-results-of-cnv-presence

This scripts allow you to filter false positives results of CNVkit program when you have not health tissue samples.
To run it you must have only python3, matplotlib, numpy.
This methos based on Bayesian model and showed good theoretical results, but this methos is not tested on real people now.

To start you should make some preparations:

Firstly you should have coverage files from "bedtools coverage" for all sample files. 
All files shuold have normalized coverage, otherwise you should create file "files_length.txt" with count of reads for every bedfile for every sample.
After that you must change file "run_filter":
1) you should have folder with  cnv probablity for every gen of cancer kind and this folder should have name "cnv_probabilities".
2) also all existed coverage files must be in folder "bed_files"
3) for you tumor sample you should specify file "inf" when cancer kind and coverage file name should be written.
File should contain:
<cancer_kind>
<coverage bed file name>
  
Note, cancer_kind must be named as cnv probabilty file in cnv_probobalites folder.

4) You can run scripts without cns file of CNVkit, but it is strongly recommended to use cns file to more accurate results. 


So, you can start filering - just run "start_filter" script. 




Requirements:
python3
matplotlib
numpy

