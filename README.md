# Filtering-false-positive-results-of-cnv-presence

These scripts allow you to filter false positive results of CNVkit tool when you have not healthy tissue samples.
CNVkit - http://cnvkit.readthedocs.io/en/stable/ . Article about CNVkit - https://www.ncbi.nlm.nih.gov/pubmed/27100738 .
To run filtering scripts you shuold have only python3, matplotlib and numpy. 
This method based on Bayesian model and showed good theoretical results, but it is not tested on real people now.

To start you should make some preparations:

<b>Firstly</b>, you should have coverage files for every reference file and sample file. Coverage files must have coverage for every bin in target bed file. All files shuold have normalized coverage, otherwise you should create file "files_length.txt" with count of reads for every bedfile. Also all coverage files must be in folder "bed_files". You can use http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html .

<b>Secondly</b>, there must be folder with name "cnv_probabilities" with "probability files", these files must contain cnv probability for every target gen for all cancer kind you are interested in.


<b>Thirdly</b>, for your tumor sample you should specify file "inf" when cancer kind and coverage file name must be written. 
File should contain:


<i>cancer_kind</i>


<i>coverage bed file name </i>
  
  
Note, cancer_kind must be named as cnv probabilty file in cnv_probobalites folder.

<b>Fourthly</b>, you can run scripts without cns file after CNVkit, but it is strongly recommended to use cns file for more accurate results. 

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


Now you can start filering - just run "start_filter" script. In output file  "result" you wiil see probability of cnv 
for every region in cns file.


<b>Requirements</b>:

python3

matplotlib


numpy


Ubuntu 16.04, but it was not tested on anything else.
