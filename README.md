# Filtering-false-positive-results-of-cnv-presence

These scripts allow you to filter false positive results of CNVkit tool when you have not health tissue samples.
To run it you must have only python3, matplotlib and numpy. 
This method based on Bayesian model and showed good theoretical results, but it is not tested on real people now.

To start you should make some preparations:

Firstly, you should have coverage files for every reference file and sample file. Coverage files must have coverage for every bin in target bed file. All files shuold have normalized coverage, otherwise you should create file "files_length.txt" with count of reads for every bedfile bedfile. Also all coverage files must be in folder "bed_files"

Secondly, there must be folder with name "cnv_probabilities" with "probability files", these files must contain cnv probability for every target gen for all cancer kind you are interested in.


Thirdly, for your tumor sample you should specify file "inf" when cancer kind and coverage file name must be written. 
File should contain:
<cancer_kind>
<coverage bed file name>
  
Note, cancer_kind must be named as cnv probabilty file in cnv_probobalites folder.

Fourthly, you can run scripts without cns file after CNVkit, but it is strongly recommended to use cns file for more accurate results. 

So, folder tree must be like that:
<pre>
project folder
├── bed_files
│   └── coverage_file1
|   └── coverage_file2 
│       
└── cnv_probabilities
|   ├── ovary
|   └── stomach
|
└── src
└── inf
└── files_length.txt
└── cns_file
└── start_filter
</pre>


So, you can start filering - just run "start_filter" script.


Requirements:


python3


matplotlib


numpy


Ubuntu 16.04, but it was not tested on something other.
