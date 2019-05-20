# mycut
A Primer Parser for NGS Sequence Data


mycut is a python script used for parsing primers from NGS data.  Inputted into the program is one FASTA file as well as the FORWARD and REVERSE primer sequences.  mycut will scan the FASTA reads for the primer sequences and, if found, will parse them out. The program takes 7 mandatory arguments:

-in_file [FILE]: Specifies Inputted FASTA File

-out_file [FILE]: Specifies Trimmed FASTA Reads Destination File

-unk_file [FILE]: Specifies Untrimmed FASTA Reads Destination File

-n_mismatch [INT]: Specifies the Allowed Number of Mismatches in Each Primer Region

-min_length [INT]: Specifies the Minimum Length Requirement for a Trimmed Read

-forward [SEQUENCE]: Specifies the Forward Primer Sequence

-reverse [SEQUENCE]: Specifies the Reverse Primer Sequence


mycut will then output three files: 

  1. one file with the trimmed and processed reads
  2. one with untrimmed reads that wern't (for whatever reason) processed
  3. a logfile
  
 


EXAMPLES:

python3 mycut_j-blamer.py –in_file merged.fasta  –out_file output.fasta –unk_file unk.fasta –n_mismatch 3 –min_len 150 -forward GTGCCAGCMGCCGCGGTAA  -reverse ACAGCCATGCANCACCT

python3 mycut_j-blamer.py –in_file  –out_file trimmed.fa –unk_file unkown.fa –n_mismatch 1 –min_len 110 -forward GGGACAATNCCACTTANAC -reverse GGCCGATTMGATTCGA
