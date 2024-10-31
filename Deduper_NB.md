Kenny's Notebook for Deduper Assignment (BI624) [on talapas...]

Assignment: https://github.com/klai22/Deduper-klai22

Data / Notebook/ WD Location: /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22

Creating interactive sessions: srun --account=bgmp -p bgmp -N 1 -c 4 --mem=100G -t 6:00:00 --pty bash

SLURM TRACKER Template: 
| Sample | JobID | Run Time(mm:ss) | CPU Usage (%)| Exit Status |
|---|---|---|---|---|
|6_2D|15903418|1:13.24|671|0|
____________________________________
10/11/24&10/14/24
____________________________________
# Part 1 – Pseudocode 

### Pseudocode 
* High-lvl Fxn (all steps below): **Dedup_SAM()**

1. Specify input arguments (input SAM file [sorted, mapped reads only], output SAM file) 
2. Create empty **set** to store  *unique* read "keys" (UMI,POS,CHROM,Strand)
3. Open input SAM file (reading)
4. Open output SAM file (writing) 
5. Iterate through EACH LINE of SAM file (to avoid loading everything into mem.)
        - '~_~' indicates a temp_variable has been created to hold value for current line 

    a. If ~line~ = header line (starts w/ '@'), write to ouput SAM file. 
    
    b. If ~line~ = read line:

        * Extracting ingredients for read's unique key 
        I. Extract ~UMI~ from QNAME column. (set # = to a temp. variable) 
            - QNAME = col1
            - UMI = everything after last ":" of QNAME string
            - High-lvl Fxn: **extract_umi()**
        II. Extract CIGAR string
            - CIGAR = col6 
        III. Calc. (true) ~POS~ so that it adjusts for soft-clipping (according to CIGAR string) 
            - Extract POS 
                - POS = col4 (set # = to a temp. variable)
            - IF CIGAR string (from ii) has an 'S' near the start (ex: 2S12M), SUBTRACT the # in front of the 'S' (2) from initial 'POS' #. 
        IV. Extract ~STRAND~ (rev_comp = TRUE [-]or FALSE[+])
            - extract ~FLAG~ = col2
            - High-lvl Fxn: **extract_strand()**
                - convert ~FLAG~ value --> convert to integer (binary)
                - Determine if read is on reverse strand: 
                    - isolate 4th bit (16 / 0x10) of converted~FLAG~ 
                    - if there is 'rev_comp', set ~strand~ = "-"
                    - else, set ~strand~ = "+"
                - set output of **extract_strand()** = ~STRAND~ 
        V. Extract ~CHROM~ 
            - RNAME = col3 (RNAME)

        * Populating set with unique keys only (& writing only unique reads to output SAM file)
        VI. Check if UMI is valid 
            -  If UMI is in list of known barcodes (STL96.txt), proceed to next steps 
            - If UMI is NOT in list of known barcodes, discard read (break & move onto next line, does not get written into output)
            - High-lvl Fxn: **validate_umi()**
        VII. Create a 'key' using ~UMI~,~POS~, ~STRAND~, & ~CHROM~
            - key can be a tuple (UMI, POS, STRAND, CHROM)
            - High-lvl Fxn: **populate_key()**
        VIII. Conditionally populate set w/ keys 
            - If the key is NOT in the set....
                1. add key to set!
                2. write read line (~line~) to output SAM file!
            - If the key is ALREADY in the set...
                - discard read (break, move onto next line, DO NOT add to set, DO NOT write into output)
6. Close input & output SAM files. 
7. exit 



### High-Level Functions 
**extract_umi**

def extract_umi(qname: str) -> str:
    
    '''
    Extracts UMI seq. from QNAME string. Appends everrything after the last ':' in QNAME. 

    REQUIRES NA (no extra fxns from bioinfo.py module needed)

    Input:
	    QNAME (string) - string extracted from col1 of given read line 
    Output:
	    UMI (string) - sequence unique to each read (unless duplicated)
    Returns extracted UMI as string
    '''

**validate_umi()**

def validate_umi(UMI: str, known_umis: set) -> bool:
    
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

**extract_strand()**

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

**populate_key()**

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

**Dedup_SAM()**

def Dedup_SAM(input: str, output: str,known_umis: set): 
    
    '''
    MAIN FUNCTION that processes SAM file (one line at a time) --> identifies PCR duplicates --> write only non-duplicate reads to output SAM file (filtering out duplicates). 

    REQUIRES NA (no extra fxns from bioinfo.py module needed)

    Input:
	    input - input SAM file. Assumes it is (1) already sorted & (2) includes ONLY MAPPED READS. 
        output - output SAM file. Same as input SAM file, but all read-lines for PCR duplicates are filtered out. 
        known_umis - a set created from STL96.txt (list of known UMIs)
    Output:
	    NA
    This function has no return value, but instead writes contents into the output SAM file. 
    '''


#### Test Examples 
**extract_umi**
```
QNAME = "NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC"
extract_umi(QNAME)
```
* Output: "CTGTTCAC" 

**validate_umi()**
```
UMI = "CTGTTCAC"
Known_UMIs={"ACGACTTG","TTCGCCTA ","ATCCATGG"}
```
* Output: False 

**extract_strand()**
Ex.1 
```
FLAG = 16
extract_strand(FLAG)
```
* Output: "-" (bc there was a 1 @ the rev_comp position)
Ex.2 
```
FLAG = 0
extract_strand(FLAG)
```
* Output: "+" (bc there was a 0 @ the rev_comp position)

**populate_key()**
```
populate_key("TGAGTGAG",52159545, "+", 2)
```
* Output: ("TGAGTGAG",52159545, "+", 2)

**Dedup_SAM()**
```
Dedup_SAM("test.sam","test_output.sam","STL96.txt")
```
* Output: NA. However, "test_output.sam" should show results indicating that PCR duplicates were correctly filtered out by function. 


SELF NOTES
* don't forget to adjust for soft clipping --> if cigar string has an s at the start, need to subtract the number in front of the S from the pos # (in sam file $ POS), 
* in the pseudocode notes, note that when you collect information for the set rmbr to (1) adjust pos for soft clipping bfore putting into set & (2) include chromosome #, & strand specificity into the details that need to be shared to ID PCr duplicates (bfore comparing UMIs)

* I created a new file &  moved all of the psuedocode into psuedocode.md 


____________________________________
10/15/24
____________________________________
# Part 1 – Pseudocode 

### Make test sam files (input & output) 
"Write examples:
* Include a properly formatted input sam file (an example sam file is included in the repo for your reference)
* Include a properly formatted expected output sam file
* Cover several different cases of things that are and are not PCR duplicates
* It may be helpful to create a "unit test folder", that contains it's own readme.md describing the test cases in your unit test samfile"


* Layout Below 

FILTERED OUT 
* A & B 
    * same UMI, CHR, POS, STRAND
* K & L & M 
    * same UMI, CHR, POS, STRAND [testing to see if it can filter out more than 2]

* N & P 
    * same UMI, CHR, STRAND
    * diff POS, BUT CIGAR STRING HAS a "S" IN IT. DO THE MATH so THAT the POS ends up summing to actually be the same AFTER CONSIDERING SOFT CLIPPING 
* Q & S 
    * same UMI, CHR, POS, STRAND
    * UMI NOT in known UMI list. 

KEPT 
* C & D
    * same CHR, POS, STRAND
    * diff UMI 
* E & F 
    * same UMI, POS, STRAND
    * diff CHR
* G & H 
    * same UMI, CHR, STRAND
    * diff POS
* I & J 
    * same UMI, CHR, POS
    * diff STRAND


* Creating Folder to hold unit test files 
```
mkdir test_files
nano test_input.sam
nano test_output.sam
nano readme.md 
```

____________________________________
10/28/24 -  P.3 -Deduper Code 
____________________________________
* fixed pseudocode.md to account for negative strands in POS calculations (when considering CIGAR string). 

* generating .py script to hold primary deduper code: 
--> /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/Lai_deduper.py



[For the documenting memory, you will need to create a bash script that runs a job for the Lai_deduper.py script]
(be prepared to report the amount of memory (GB) your script required and the time it took to run!)
* will need to report # of header lines, unique reads, wrong UMIs, & duplicated removed. 
* will also need to report the # of reads per chrom. in tsv format. 
* Will run on "/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam"



* Updating my .test files to account for (-) CIGAR STRING EXAMPLES 

* Edited the CIGAR string below to be more complex, changed FLAG to be (+) 
# N & P (same UMI, CHR, STRAND, --> different POS, but CIGAR with soft clipping makes POS the same) AATTCCGG (+)
N:X:X:1:1:1:1:AATTCCGG  163  chr1    100 60  4S10M   =   400 0 GCAA...GCAA IJKL789
P:X:X:1:1:1:1:AATTCCGG  163  chr1    99  60  3S5M2D3M6S  =   400 0   GCAA...GCAA IJKL789
* Added pair below for CIGAR string and (-) FLAG 
# O & R (same UMI, CHR, STRAND, --> different POS, but CIGAR with soft clipping makes POS the same) CAACTGGT (-) 
O:X:X:1:1:1:1:CAACTGGT  83 chr1 230 60  2S3M1I5N1D1S    =   200 0   AAGT...GGCG ABCD123
R:X:X:1:1:1:1:CAACTGGT  83 chr1 228 60  6M7I3D3S    =   200 0   TACG...CCAT ABCD123



* Re-installing samtools because I don't remember which conda env I put it in. 
```
conda create --name bgmp_samtools 
conda activate bgmp_samtools
conda install -c bioconda samtools
samtools --version
```
* Sorting my test sam file: 
```
cd /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/test_files
samtools sort -n -o test_input_sorted.sam --output-fmt sam test_input.sam

```

ERROR: 
[W::sam_read1] Parse error at line 25
samtools sort: truncated file. Aborting
____________________________________
10/30/24 -  P.3 -Deduper Code 
____________________________________
* Fixing my test .sam file : 
    * /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/test_files/test_input.sam
        * I am changing all chromosome values in CHR field from 'chr1'-->'1' to match the headers at the top (Leslie said my chromsome names do not match )

* Reatempting samtools sort 
```
cd /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/test_files
conda activate bgmp_samtools
#sorting 
samtools sort -o test_input_sorted.sam --output-fmt sam test_input.sam
```
* It finally sorted correctly!
--> CURRENT TEST SAM FILE TO USE: /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/test_files/test_input_sorted.sam


* I wrote all high-level functions into /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/Lai_deduper.py 
```
#Code Ran 
./Lai_deduper.py -f "test_files/test_input.sam" -o "empty.sam" -u "STL96.txt"
#Output 
extract_umi() working properly
validate_umi() working properly
extract_strand() working properly 
populate_key() working properly
calc_pos() working properly
```

____________________________________
10/31/24 -  P.3 -Deduper Code 
____________________________________
* Working on Code to Load File in /home/kenlai/bgmp/bioinfo/Bi624/Assignments/Deduper-klai22/Lai_deduper.py
    * ended on validate_umi(), next steps are to code 'VII. Create a 'key' using ~UMI~,~POS~, ~STRAND~, & ~CHROM~'
    * NEXT, also need to test my code so far, I haven't been testing what I wrote so far. 








TO DO: 
* Edit script (as your issues suggested) to use sets instead of tuples (ppl were against the idea of a "key" to avoid time consumption)
    https://github.com/klai22/Deduper-klai22/issues
* Create a .sh script that will get metrics such as GB used by Lai_deduper.py & also includes steps like 'samtools sort'-ing the files before hand. 
________________________

# template stuff 
*  Createded a conda env called QAA  & installed FASTQC 
```
$ conda create --name QAA
$ conda activate QAA
$ conda install bioconda::fastqc


# SAM STRUCTRURE REMINDER: 
QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL