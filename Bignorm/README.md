This is the code for the Bignorm Illumina read normalization tool.

You'll need a C compiler with OpenMP (tested with gcc >= 4.8) and 
a development version of libz (not too old or the '-z' flag 
of Bignorm will be ignored and all output will be gzip'ed)

Just compile using 'make'.

Basic usage:
let reads_1.fq.gz and reads_2.fq.gz the input files to normalize, then
the simplest command is:

Bignorm -1 reads_1.fq.gz -2 reads_2.fq.gz

which will produce reads_1.fq.gz_keep and reads_2.fq.gz_keep.

For advanced options, try

Bignorm -?

or read the manual

