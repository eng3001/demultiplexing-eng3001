#!/bin/python
import gzip
import argparse
import re

#Set input variables
def get_args():
    parser = argparse.ArgumentParser(description="Generate a per base \
    distribution of quality scores for a file")
    parser.add_argument("-f", "--file", help="File name", type=str)
    parser.add_argument("-r", "--run", help="Run File Number", type=str)
    return parser.parse_args()

args = get_args()

#Set global variables
FILE = args.file
R_var = args.run

#All file names, dont need just yet
#index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
#R1_file = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
#R2_file = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
#R3_file = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
#R4_file = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

index_dict = {}
sum_q_array = []
rec_count = 0
#line count = 1452986940

def convert_phred(letter):
    """Converts a single character into a phred score"""
    val = ord(letter)
    val = val - 33
    return val

#Code to load dictionary with known indexes
#with open(index_file, "r") as ifh:
#    ifh.readline()
#    for line in ifh:
#        line = line.split()
#        index = line[4]
#        index_dict.setdefault(index, 1)

with gzip.open(FILE, 'rt') as zip:
    LN = 0
    for line in zip:
        LN+=1
        line = line.strip()
        if LN == 4:
            length = len(line)
            for x in range(length):
                sum_q_array.append(0)
        if LN % 4 == 0:
            rec_count+=1
            counter = 0
            for letter in line:
                #print (letter)
                sum_q_array[counter] += (convert_phred(letter))
                counter+=1

print (sum_q_array)

index = 0
while index < len(sum_q_array):
    sum_q_array[index] = (sum_q_array[index]/rec_count)
    index += 1

index_array = []
index = 1
while index <= len(sum_q_array):
    index_array.append(index)
    index+=1

import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
plt.plot(index_array,sum_q_array, "o")
plt.xlabel('Base Pair')
plt.ylabel('Mean Q-score')
plt.title('Distribution of Q-score')
plt.savefig("%s.png" % (R_var))
#plt.yscale('log')
#plt.xlim(0,10000)
#plt.show()

#print (index_array)
#print (sum_q_array)
#print (rec_count)
#print ('A = ', convert_phred('A'))
#print ('3 As =', convert_phred('A')*3)
