# Demultiplexing

The objective of this function is to demultiplex Illumina reads that are dual-indexed. The script can identify low quality score reads and index hopping. Script will produce R1 and R2 fastq files for each known index, unknow records, and index hopped records.

### Input Parameters
- -f1,--R1_file : Read 1 File name.
- -f2,--R2_file : Read 2 File name. 
- -f3,--R3_file : Read 3 File name.
- -f4,--R4_file : Read 4 File name.
- -i, --index_file : Index file name.
- -iq,--index_cutoff : Minimum average index quality score cut off

Required parameter to run ```Demultiplexing_Checking_Sequence_Quality_Score_longer_runtime.py```
- -sq,--seq_cutoff : Minimum average sequence quality score cut off
