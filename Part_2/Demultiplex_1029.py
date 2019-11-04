#!/bin/python
import gzip
import argparse
import re
import numpy as np
from itertools import islice


#-----------------------------------------------------------------------
# The objective of this function is to demultiplex Illumina reads that are
# dual-indexed. The script can identify low quality score reads and
# index hopping. Script will produce R1 and R2 fastq files for each known
# index, unknow records, and index hopped records.
#-----------------------------------------------------------------------


#Set input variables
def get_args():
    parser = argparse.ArgumentParser(description="Generate a per base \
    distribution of quality scores for a file")
    parser.add_argument("-f1", "--R1_file", help="R1 File name", type=str)
    parser.add_argument("-f2", "--R2_file", help="R2 File name", type=str)
    parser.add_argument("-f3", "--R3_file", help="R3 File name", type=str)
    parser.add_argument("-f4", "--R4_file", help="R4 File name", type=str)
    parser.add_argument("-i", "--index_file", help="index file name", type=str)
    parser.add_argument("-sq", "--seq_cutoff", help="minimum average sequence quality score cut off", type=int)
    parser.add_argument("-iq", "--index_cutoff", help="minimum average index quality score cut off", type=int)
    return parser.parse_args()

args = get_args()
#Given Files
R1_seq_file = args.R1_file
R1_index_file = args.R2_file
R2_index_file = args.R3_file
R2_seq_file = args.R4_file
index_file = args.index_file
#Minimum quality score input values
seq_qual_score = args.seq_cutoff
index_qual_score = args.index_cutoff

#  Will be used to check qscore.
def convert_phred(letter):
    '''Converts a single character into a phred score'''
    if letter != 'N':
        val = ord(letter)
        val = val - 33
    else:
        val = 0
    return val

#Used in reverse complement function
def reverse(string):
    '''Reverse a string. Used in reverse_complement'''
    str = ""
    for char in string:
        str = char + str
    return str

def reverse_complement(sequence):
    '''Finds the reverse complement of a given string'''
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    temp_string = str()
    for character in sequence:
        temp_string += complement_dict[character]
    return reverse(temp_string) #Returning the reverse of the complement

def write_record(file_name, record_array, barcode1, barcode2):
    '''Writes a record to the corresponding file'''
    #File names stored as the value of the dictionary
    with open (file_name, "a") as write_file:
        header_line = record_array[0] + "_" + barcode1 + "_" + barcode2
        write_file.write(header_line + '\n')
        write_file.write(record_array[1] + '\n')
        write_file.write(record_array[2] + '\n')
        write_file.write(record_array[3] + '\n')
        write_file.close()

def check_index_quality(index_string, seq):
    '''Checks if the index q-scores are above a specified value'''
    sum_q_array = []
    bool_val = False
    length = len(index_string)
    for x in range(length):
        sum_q_array.append(0)
    counter = 0
    for letter in index_string:
        sum_q_array[counter] += (convert_phred(letter))
        counter+=1
    numpy_array = np.array(sum_q_array)
    mean = numpy_array.mean()
    if mean >= index_qual_score: #Makes sure that the quality score is above a average of 28
        bool_val = True
    if 'N' in seq:
        bool_val = False
    return bool_val

def check_seq_quality(seq):
    '''Checks if the sequences average q-score is above a specific value'''
    sum_q_array = []
    length = len(seq)
    for x in range(length):
        sum_q_array.append(0)
    bool_val = False
    counter = 0
    for letter in seq:
        sum_q_array[counter] += (convert_phred(letter))
        counter+=1
    numpy_array = np.array(sum_q_array)
    mean = numpy_array.mean()
    if mean >= seq_qual_score: #Makes sure that the average quality score is 35 or above
        bool_val = True
    return bool_val

def valid_barcode(barcode, barcode_dict):
    '''Checks to see if the sequenced barcode is in the given barcode list'''
    bool_val = False
    if barcode in barcode_dict:
        bool_val = True
    return bool_val

def read_record(file_handle):
    '''Reads a record from a file and outputs it as a list'''
    rec_array = []
    record_object = islice(file_handle, 4) #islice returns an iterable object
    rec_list = list(record_object)
    for line in rec_list:
        line=line.strip()
        rec_array.append(line)
    return(rec_array)

#Dictionaries used to store the indexes
R1_dict = dict()
R2_dict = dict()

#Dictionary to keep track of the counts of each index
index_count_dict = dict()

#Given Files
#R1_seq_file = "/projects/bgmp/wyatte/Bi622/1294_S1_L008_R1_001.fastq.gz"
#R1_index_file = "/projects/bgmp/wyatte/Bi622/1294_S1_L008_R2_001.fastq.gz"
#R2_index_file = "/projects/bgmp/wyatte/Bi622/1294_S1_L008_R3_001.fastq.gz"
#R2_seq_file = "/projects/bgmp/wyatte/Bi622/1294_S1_L008_R4_001.fastq.gz"

# Files to output Unknown & Index Hopped reads
Unknown_R1 = "Unknown_R1.fq"
Unknown_R2 = "Unknown_R2.fq"
Hopped_R1 = "Index_Hopped_R1.fq"
Hopped_R2 = "Index_Hopped_R2.fq"

#Read index file, create dictionary with key being indexes and value being file name
with open("/projects/bgmp/wyatte/Bi622/indexes.txt") as fh:
    fh.readline()
    for line in fh:
        line = line.strip().split()
        index = line[4]
        R1_dict[index] = index + "_R1.fq"
        R2_dict[index] = index + "_R2.fq"

#Body of code
with gzip.open(R1_seq_file, 'rt') as R1_seq, gzip.open(R1_index_file, 'rt') as R1_index, gzip.open(R2_index_file, 'rt') as R2_index, gzip.open(R2_seq_file, 'rt') as R2_seq:
    line = "Hi"
    read_pairs_val = 0
    index_hopped_val = 0
    unknown_pairs_val = 0
    total_recs = 0

    while True:
        seq_R1_qlty = False #Used to check if the quality of the sequence is above standard
        seq_R2_qlty = False
        index_R1_qlty = False #Used to check if the quality of the index is above standard
        index_R2_qlty = False

        #Checks to see if the indexes are valid
        Valid_Index_R1 = False
        Valid_Index_R2 = False

        #Arrays to store a record from each file
        R1_seq_array = []
        R1_index_array = []
        R2_seq_array = []
        R2_index_array = []

        #Read in 4 lines from each file at a time
        R1_seq_array = read_record(R1_seq)
        R1_index_array = read_record(R1_index)
        R2_index_array= read_record(R2_index)
        R2_seq_array = read_record(R2_seq)

        #Checks to see if its the end of the file by checking the length of the records read in
        if (len(R1_seq_array) + len(R1_index_array) + len(R2_seq_array) + len(R2_index_array)) < 16:
            break

        # Checking to see if the quality scores for the indexes and sequences are good
        seq_R1_qlty = check_seq_quality(R1_seq_array[3])
        seq_R2_qlty = check_seq_quality(R2_seq_array[3])
        index_R1_qlty = check_index_quality(R1_index_array[3], R1_index_array[1])
        index_R2_qlty = check_index_quality(R2_index_array[3], R2_index_array[1])

        #Increase the of records to keep track of the total
        total_recs += 1

        # Checking to see if the indexes are valid
        Valid_Index_R1 = valid_barcode(R1_index_array[1], R1_dict)
        Valid_Index_R2 = valid_barcode(reverse_complement(R2_index_array[1]), R2_dict)

        #Chunk of code that filters through the records and prints it to the right file
        if (seq_R1_qlty & seq_R2_qlty & index_R1_qlty & index_R2_qlty & Valid_Index_R1 & Valid_Index_R2):
            if R1_index_array[1] == reverse_complement(R2_index_array[1]):
                read_pairs_val += 1
                write_record(R1_dict[R1_index_array[1]], R1_seq_array, R1_index_array[1], R2_index_array[1])
                write_record(R2_dict[reverse_complement(R2_index_array[1])], R2_seq_array, R2_index_array[1], R1_index_array[1])
                if R1_index_array[1] in index_count_dict.keys():
                    index_count_dict[R1_index_array[1]] += 1
                else:
                    index_count_dict.setdefault(R1_index_array[1],1)
            else:
                index_hopped_val += 1
                write_record(Hopped_R1, R1_seq_array, R1_index_array[1], R2_index_array[1])
                write_record(Hopped_R2, R2_seq_array, R2_index_array[1], R1_index_array[1])
        else:
            unknown_pairs_val += 1
            write_record(Unknown_R1, R1_seq_array, R1_index_array[1], R2_index_array[1])
            write_record(Unknown_R2, R2_seq_array, R2_index_array[1], R1_index_array[1])

    #print summary report
    print("Total number of records processed: " + str(total_recs))
    print("Number of read pairs with properly matched indexes: " + str(read_pairs_val))
    print("Number of read pairs with index hopping observed: " + str(index_hopped_val))
    print("Percent of records where index hopping was observed: " + str((index_hopped_val/total_recs)*100) + "%")
    print("Number of read pairs with unknown indexes: " + str(unknown_pairs_val))
    print("Percent of records with unknown indexes: " + str((unknown_pairs_val/total_recs)*100) + "%")
    print("-------------------------------------------------------------")
    print("INDEX | Count of Index Reads | Percentage of Index Reads")
    for key in index_count_dict.keys():
        print(str(key) + " | Count: " + str(index_count_dict[key]) + " | " + str((index_count_dict[key]/total_recs)*100) + "%")

    #Closing Files
    R1_seq.close()
    R1_index.close()
    R2_seq.close()
    R2_index.close()
