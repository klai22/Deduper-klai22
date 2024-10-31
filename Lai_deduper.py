#!/usr/bin/env python

#./sam_parse.py -sam_file ___

#./Lai_deduper.py -f "test_files/test_input_sorted.sam" -o1 "filtered_sorted.sam" -o2 "PCR_dups.tsv" -o3 "removed_reads.tsv" -u "STL96.txt"

# MANAGING INPUTS 
#setting get_args 
import argparse
#defining input args 
def get_args():
    parser = argparse.ArgumentParser(description="Deduper script to remove PCR duplicates from a SAM file. This code assumes (1) a sorted sam file & (2) single-end reads only/ not paired-end")
    parser.add_argument("-f", required=True, help="Absolute path to INPUT sorted .sam (pre-dedup.)") 
    parser.add_argument("-o1", required=True, help="Absolute path to OUTPUT sorted .sam (post-dedup.)")
    parser.add_argument("-o2", required=True, help="Absolute path to OUTPUT aligned reads that were identified as PCR Duplicates")
    parser.add_argument("-o3", required=True, help="Absolute path to OUTPUT aligned reads that were removed for reasons besides PCR duplication (such as invalid UMIs)")    
    parser.add_argument("-u", required=True, help="Absolute path to .txt of known UMIs")
    return parser.parse_args()

#setting input args--> variables 
args=get_args()
input_sam_file=args.f
output_sam_file=args.o1
PCR_duplicates_file=args.o2
removed_reads_file=args.o3
known_umis_file=args.u

#Importing modules
# import gzip
# import matplotlib.pyplot as plt
# import numpy as np
import re # a module that provides regex matching operations https://docs.python.org/3/library/re.html#re.search

# HIGH-LVL FXNs 
# **extract_umi()**
def extract_umi(qname: str) -> str:
    '''
    Extracts UMI seq. from QNAME string. Appends everything after the last ':' in QNAME. 
    REQUIRES NA (no extra fxns from bioinfo.py module needed)
    Input:
        QNAME (string) - string extracted from col1 of given read line 
    Output:
        UMI (string) - sequence unique to each read (unless duplicated)
    Returns extracted UMI as string
    '''
    #Split qname by ':', isolate last element (the UMI)
    UMI = qname.split(':')[-1]
    return UMI

# **validate_umi()**
def validate_umi(UMI:str, known_umis:set)-> bool:
    '''
Checks if UMI is also present in list (or set) of known UMIs. Error correction could be encorperated into this function w/ hamming distance thresholds if desired. 

REQUIRES NA (no extra fxns from bioinfo.py module needed)

Input:
    UMI (string) - output of extract_umi()
    known_umis - a set created from STL96.txt (list of known UMIs)
Output:
    validation result (bool) - T or F (whether UMI is known list or not)
Returns True if UMI is in known_UMI list, False if no matches detected. 
    '''
    return UMI in known_umis

# **extract_strand()**
def extract_strand(FLAG: int) -> str:
    '''
Extracts strandedness given ~FLAG~ (isolated col2 of read line).

REQUIRES NA (no extra fxns from bioinfo.py module needed)

Input:
    FLAG (string) - integer, isolated 2nd field of read line
Output:
    strandedness(str) - "+" or "-" depending on presence of rev_comp in flag 
Returns "+" if no rev_comp present in flag. Returns "-" if rev_comp detected in flag. 
    '''
    #Check if the 16th bit is set 
    if FLAG & 0x10: #0x10 = 16 in 'hexicimal', refers to 16th bit 
        return "-" #if T, then it is a RVS strand 
    else:
        return "+" #if F, then it is on FWD strand 

# **populate_key()**
def populate_key(UMI: str, POS: int,STRAND: str, CHROM: int)-> tuple:
    '''
Compiles the 4 main variables of a given read line that we compare to determine if read is a PCR duplicate. 

REQUIRES NA (no extra fxns from bioinfo.py module needed)

Input:
    UMI (string) - output of extract_umi()
    POS (int) - starting pos. of read (factoring in calc. for soft-clipping!)
    STRAND (str) - output of extract_strand()
    CHROM (int) - RNAME value for given read line 
Output:
    key (tuple) - a key to be conditionally encorperated into the initalized set (depending on whether or not duplicate exists)
Returns a tuple of UMI, POS, STRAND, & CHROM for given read line. This tuple should be unique to the read if it is NOT a PCR duplicate. 
'''
    return (UMI,POS,STRAND,CHROM)

# **calc_pos()**
def calc_pos(POS: int,STRAND: str, CIGAR: str)-> int:
    '''
Calculates the "true POS" (col4) considering CIGAR string & strandedness. 
    * IF [+]
    IF CIGAR string has an 'S' near the start (ex: 2S12M), SUBTRACT the # in front of the 'S' (2) from initial 'POS' #.
    * IF [-] 
    ADD M,D,N,&(right-most ONLY)S to initial 'POS'#. 
    
REQUIRES NA (no extra fxns from bioinfo.py module needed)

Input:
    POS (int) - starting pos. of read (factoring in calc. for soft-clipping!)
    STRAND (str) - output of extract_strand()
    CIGAR (str) - CIGAR string (ex:5S10M2I3D) 
Output:
    Updated_POS (int) - True position # (shared b/w true PCR duplicates) calc. by considering factors like soft-clipping in the CIGAR string 
Returns an integer of what the updated POS should be. THIS will be stored in the tuples/keys. 
    '''
    #Initalize the Updated POS # 
    new_pos=POS
    
    #If strand is + 
    if STRAND == "+":
        # Find left-most 'S' 
            # ('re' module allows us to us regex in pyt., re.match(pattern, string) ONLY CHECKS FOR A MATCH @ BEGINNING OF STRING --> store that # to S_count) 
            # ' r"(..)" ' indicates it is a raw string --> \ treated no longer an escape chr. 
        S_count=re.match(r"([0-9]+)S",CIGAR)
            # If there is S @ beginning of string, SUBTRACT from POS
        if S_count:
            #The actual # in front of S needs to be extracted from the created S_count obj. 
            S_value = int(S_count.group(1))
            new_pos-=S_value
    #If strand is - 
    elif STRAND == "-":
        #Find all relevant elements of CIGAR 
            #'M'atch/mismatch, 'D'eletions, & Skipped ('N')
                #findall -> extracts all #s preceding letters M, D, & N 
        MDN_count = re.findall(r"([0-9]+)([M,D,N])", CIGAR)
            #Right-most 'S' , re.search similar to re.match, recall '$' indicates end of (right-most occurance) the string in regex 
        S_count = re.search(r"([0-9]+)S$",CIGAR)
        
        #Sum up all counts per letter type (M,D,N) from MDN_count & add to new_pos (we don't use if here bc there should always be atleast 1M )
        for count,letter_type in MDN_count:
            new_pos += int(count)
        #Add on right-most S count to new_pos 
        if S_count:
            new_pos += int(S_count.group(1))
    return new_pos


# DEDUP CODE 

# 1. Specify input arguments (input SAM file [sorted, mapped reads only], output SAM file) 
    #Completed Earlier 
# 2. Create empty **set** to store  *unique* read "keys" (UMI,POS,CHROM,Strand)
#open the known umis file 
with open(known_umis_file,"r") as file:
    # strip /n for every line in file, turn all into a set 
    known_umis = set(line.strip() for line in file)
# 3. Open input SAM file (reading) & # 4. Open output SAM file (writing) 
with open (input_sam_file,"r") as input_sam, open(output_sam_file, "w") as output_sam, open(PCR_duplicates_file, "w") as PCR_dups, open(removed_reads_file, "W") as removed_reads: 
    # 5. Iterate through EACH LINE of SAM file (to avoid loading everything into mem.)
        #- '~_~' indicates a temp_variable has been created to hold value for current line 
    for line in input_sam:
        #a. If ~line~ = header line (starts w/ '@'), write to ouput SAM file. 
        if line.startswith("@"):
            output_sam.write(line)
            continue #skips to next line, going back to start of for loop 
        #b. If ~line~ = read line:
    #Extracting ingredients for read's unique key 
        #breaking up line (fields) into multiple accessible components - tab-delimited 
        fields = line.strip().split("/t")
        QNAME = fields[0]
        FLAG = int(fields[1])
        CHR = int(fields[2])
        POS = int(fields[3])
        CIGAR = fields[5]
        #I. Extract ~UMI~ from QNAME column. (set # = to a temp. variable)
        UMI=extract_umi(QNAME)
        #II. Extract ~STRAND~ (rev_comp = TRUE [-]or FALSE[+])
        STRAND=extract_strand(FLAG)
        #III. Extract CIGAR string & IV. Calc. (true) ~POS~ so that it adjusts for soft-clipping & other elements (according to CIGAR string) 
        NEW_POS = calc_pos(POS, STRAND, CIGAR)
        #V. Extract ~CHROM~ 
            #Completed Earlier
    #Populating set with unique keys only (& writing only unique reads to output SAM file)
        #VI. Check if UMI is valid
        if validate_umi(UMI, known_umis)==False:
            removed_reads.write(line) #write invalid umis into removed_reads file (read that were removed that weren't PCR duplicates)
            continue #move to next read line, going back to for loop again 
        # DONT NEED 'elif validate_umi(UMI, known_umis)==True:', script will know to proceed directly to next steps 


#         VII. Create a 'key' using ~UMI~,~POS~, ~STRAND~, & ~CHROM~
#             - key can be a tuple (UMI, POS, STRAND, CHROM)
#             - High-lvl Fxn: **populate_key()**
#         VIII. Conditionally populate set w/ keys 
#             - If the key is NOT in the set....
#                 1. add key to set!
#                 2. write read line (~line~) to output SAM file!
#             - If the key is ALREADY in the set...
#                 - discard read (break, move onto next line, DO NOT add to set, DO NOT write into output)
# 6. Close input & output SAM files. 
# 7. exit 





# UNIT TESTS 
assert extract_umi("NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC") == "CTGTTCAC", f"extract_umi() working properly"
print("extract_umi() working properly")
test_UMI = "CTGTTCAC"
test_Known_UMIs={"ACGACTTG","TTCGCCTA ","ATCCATGG"}
assert validate_umi(test_UMI, test_Known_UMIs) == False, f"validate_umi() fao;ed"
print("validate_umi() working properly")
assert extract_strand(16) =="-", f"extract_strand() failed"
assert extract_strand(0) =="+", f"extract_strand() failed"
assert extract_strand(147) =="-", f"extract_strand() failed"
assert extract_strand(99) =="+", f"extract_strand() failed"
assert extract_strand(83) =="-", f"extract_strand() failed"
assert extract_strand(163) =="+", f"extract_strand() failed"
print("extract_strand() working properly ")
unit_test_A=populate_key("TGAGTGAG",52159545, "+", 2)
assert unit_test_A == ("TGAGTGAG",52159545, "+", 2), f"populate_key() failed"
assert isinstance(unit_test_A, tuple) == True, f"populate_key() failed"
print("populate_key() working properly")
assert calc_pos(100,"+","4S10M") == 96,f"calc_pos() failed" #N 100-3
assert calc_pos(99,"+","3S5M2D3M6S") == 96,f"calc_pos() failed" #P 99-3
assert calc_pos(230,"-","2S3M1I5N1D1S") == 240,f"calc_pos() failed" #O 230+10
assert calc_pos(228,"-","6M7I3D3S") == 240,f"calc_pos() failed" #R 228+12
assert calc_pos(3,"-","5S10M2I3D")==16, f"calc_pos() failed" #3+13
print("calc_pos() working properly")