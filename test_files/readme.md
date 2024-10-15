The information below describes the read lines in the test .sam files (test_input.sam & test_output.sam) and the conditions each pair/group were attempting to test in terms of succesfully deduplication. 

The Read Name (each assigned a letter) is indicated as the FIRST CHARACTER in the QNAME of each read line in said test .sam files. **For the columns UMI, CHR, POS, & STRAND, an X is marked if said category is annotated as DIFFERENT between "test-pairs" (pairings indicated in a column 2).** If no X is drawn in any of these categories, they are all the same between test pair (reads). 


| Read Name | Test Pair | UMI | CHR| POS | STRAND | Condition Tested | Correct Result| Known UMI (if assigned)| 
|---|---|---|---|---|---|---|---|---|
|A|B|||||Correctly filtering out PCR duplicates that share all 4 categories|Filtered Out|AACGCCAT|
|B |A |||||Correctly filtering out PCR duplicates that share all 4 categories|Filtered Out|AACGCCAT|
|K |L,M |||||Correctly filtering our >2 PCR duplicates that share all 4 categories|Filtered Out|AAGGTACG|
|L |K,L |||||Correctly filtering our >2 PCR duplicates that share all 4 categories|Filtered Out|AAGGTACG|
|M |K,M |||||Correctly filtering our >2 PCR duplicates that share all 4 categories|Filtered Out|AAGGTACG|
|N |P |||X||Technically  PCR duplicates, but have different POS. Should be classified as duplicated once soft-clipping is accounted for. |Filtered Out|AATTCCGG|
|P |N |||X||Technically  PCR duplicates, but have different POS. Should be classified as duplicated once soft-clipping is accounted for. |Filtered Out|AATTCCGG|
|Q |S |||||Succesfully filter out reads w/ invalid UMIs (UMI not detected in known UMI list)|Filtered Out|NA|
|S |Q |||||Succesfully filter out reads w/ invalid UMIs (UMI not detected in known UMI list)|Filtered Out|NA|
|C |D |X||||Correctly identifies non-PCR duplicates (different UMIs)|Kept|ACACAGAG|
|D |C |X||||Correctly identifies non-PCR duplicates (different UMIs)|Kept|TTCGTTCG|
|E |F ||X|||Correctly identifies non-PCR duplicates (different CHROMs)|Kept|ACACTCAG|
|F |E ||X|||Correctly identifies non-PCR duplicates (different CHROMs)|Kept|ACACTCAG|
|G |H |||X||Correctly identifies non-PCR duplicates (different POSs). CIGAR strings are the same|Kept|TTCGCCTA|
|H |G |||X||Correctly identifies non-PCR duplicates (different POSs). CIGAR strings are the same|Kept|TTCGCCTA|
|I |J ||||X|Correctly identifies non-PCR duplicates (different STRANDs). FLAG of 83 codes for a (-) strand. |Kept|GATCCTAG|
|J |I ||||X|Correctly identifies non-PCR duplicates (different STRANDs). FLAG of 163 codes for a (+) strand. |Kept|GATCCTAG|

* For example if Read Name is C, read line's QNAME appears as "**C**:X:X:1:1:1:1:ACACAGAG" in test .sam file. 



