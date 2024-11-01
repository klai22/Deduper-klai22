Kenneth Lai 
BI624: Deduper Part 1 - Pseudocode 

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
        II. Extract ~STRAND~ (rev_comp = TRUE [-]or FALSE[+])
            - extract ~FLAG~ = col2
            - High-lvl Fxn: **extract_strand()**
                - convert ~FLAG~ value --> convert to integer (binary)
                - Determine if read is on reverse strand: 
                    - isolate 4th bit (16 / 0x10) of converted~FLAG~ 
                    - if there is 'rev_comp', set ~strand~ = "-"
                    - else, set ~strand~ = "+"
                - set output of **extract_strand()** = ~STRAND~ 
        III. Extract CIGAR string
            - CIGAR = col6 
        IV. Calc. (true) ~POS~ so that it adjusts for soft-clipping (according to CIGAR string) 
            - High-lvl Fxn: **calc_pos()**
            - Extract POS 
                - POS = col4 (set # = to a temp. variable)
            - IF [+] 
                - IF CIGAR string (from ii) has an 'S' near the start (ex: 2S12M), SUBTRACT the # in front of the 'S' (2) from initial 'POS' #. 
            - IF [-] 
                - ADD M,D,N,&(right-most ONLY)S to initial 'POS'#. 
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

**calc_pos()**

def calc_pos(POS: int,STRAND: str, CIGAR: str)-> int: 
    
    '''
    Calculates the "true POS" considering CIGAR string & strandedness. 

    REQUIRES NA (no extra fxns from bioinfo.py module needed)

    Input:
        POS (int) - starting pos. of read (factoring in calc. for soft-clipping!)
        STRAND (str) - output of extract_strand()
        CIGAR (str) - CIGAR string (ex:5S10M2I3D) 
    Output:
	    Updated_POS (int) - True position # (shared b/w true PCR duplicates) calc. by considering factors like soft-clipping in the CIGAR string 
    Returns an integer of what the updated POS should be. THIS will be stored in the tuples/keys. 
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

**calc_pos()**
```
calc_pos(3,"-","5S10M2I3D")
```
* Output: 16 (bc 3+13(10M+3D)=16)