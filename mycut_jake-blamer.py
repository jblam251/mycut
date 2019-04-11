import os
import sys
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

## Help Screen
if sys.argv[1] == "--h":
	print("\n\n-----Thanks for Using mycut!-----\n\n")
	print("\nSYNOPSIS\n")
	print("\tpython3 mycut_j-blamer.py [ –in_file –out_file –unk_file –n_mismatch –min_len -forward -reverse ]\n")
	print("\nDESCRIPTION\n")
	print("\tmycut is a python script used for parsing primers from NGS data.  Inputted into the program is one FASTA file as well as the FORWARD and REVERSE primer sequences.  mycut will scan the FASTA reads for the primer sequences and, if found, will parse them out. mycut takes 7 mandatory arguments:\n\n\t -in_file [FILE]: Specifies Inputted FASTA File\n\t -out_file [FILE]: Specifies Trimmed FASTA Reads Destination File\n\t -unk_file [FILE]: Specifies Untrimmed FASTA Reads Destination File\n\t -n_mismatch [INT]: Specifies the Allowed Number of Mismatches in Each Primer Region\n\t -min_length [INT]: Specifies the Minimum Length Requirement for a Trimmed Read\n\t -forward [SEQUENCE]: Specifies the Forward Primer Sequence\n\t -reverse [SEQUENCE]: Specifies the Reverse Primer Sequence\n\n")
	print("\nOUTPUT\n")
	print("\tmycut will output three files: one file with the trimmed and processed reads, one with untrimmed reads that wern't, for whatever reason, processed, and a logfile with basic statistics and explinations for why unprocessed reads wern't processed\n") 
	print("\nEXAMPLES\n")
	print("\tpython3 mycut_j-blamer.py –in_file merged.fasta  –out_file output.fasta –unk_file unk.fasta –n_mismatch 3 –min_len 150 -forward GTGCCAGCMGCCGCGGTAA  -reverse ACAGCCATGCANCACCT\n")
	print("\tpython3 mycut_j-blamer.py –in_file  –out_file trimmed.fa –unk_file unkown.fa –n_mismatch 1 –min_len 110 -forward GGGACAATNCCACTTANAC -reverse GGCCGATTMGATTCGA\n")
	

## Set Arguments as Variables
elif sys.argv[1] == "-in_file":
	input_file = sys.argv[2]
	
	if sys.argv[3] == "-out_file":
		output_file = sys.argv[4]
	if sys.argv[5] == "-unk_file":
		unk_file = sys.argv[6]
	if sys.argv[7] == "-n_mismatch":
		threshold = sys.argv[8]
	if sys.argv[9] == "-min_len":
		min_length = sys.argv[10]
	if sys.argv[11] == "-forward":
		primer = sys.argv[12]
	if sys.argv[13] == "-reverse":
		r_primer = sys.argv[14]

	
	primer = Seq(primer, generic_dna)
	r_primer = Seq(r_primer, generic_dna)
	RCr_primer = r_primer.reverse_complement()

	## Set Variables and Read Fasta Lines into Dictionary for Easier Access 
	fa_open = open(input_file)
	as_list = fa_open.readlines()
	ln_count = 0
	trim_count = 0
	total_trimmed_length = 0
	total_untrimmed_length = 0
	fa_dict = {}

	## For Dictionary: Fasta Headers as Keys, Sequences as Values
	for line in as_list:
		ln_count = ln_count + 1
		if ln_count % 2 == 0:
			id = line.rstrip()
		elif ln_count % 2 != 0:
			seq = line.rstrip()
		fa_dict[seq] = id
	
	## Store Keys and Values in Variables for Easier Access
	value_list = fa_dict.values()
	key_list = fa_dict.keys()
	
	## Search for Primers, Iterate through Each Read (key)
	for key in key_list:
		primer_length = len(primer) 
		Rprimer_length = len(RCr_primer)
		read_length = len(fa_dict[key]) 
		max_itrs = read_length - Rprimer_length
		## Objects fwd and rvs Mark if Forward/Reverse Primer Detected
		fwd = 0
		rvs = 0

		for i in range(0,max_itrs+1):
			score = 0
			r_score = 0
			read = fa_dict[key]
			if rvs == 1:
				break
			
			## Allow Mismatches and Handle Ambiguities for Forward
			for j in range(0,primer_length):
				if fwd == 1:
					break
				if primer[j] == "N":
					score = score + 1
				## In BioPython K is Reverse Complement of M
				if primer[j] == "K":
					if read[i+j] == "A":
						score = score + 1	
					if read[i+j] == "C":
						score = score + 1
				## Score -1 for a Mismatch
				if read[i+j] != primer[j]:
					score = score - 1
				## Anytime Score Falls Below Threshold, Break
				if score < -(int(threshold)):
					break
				## Keep if Read Remains Above Threshold AND Spans Full Primer LEngth
				if j == primer_length - 1:
					fwd = 1
					## Cut Forward Primer
					cut_read = read[i+len(primer):]
					break
				## Break after First 50 Bases
				if i >= 50:
					break
	
			## Reverse is only searched if Foward (fwd) was Detected
			if fwd == 1:
				## For Searching Reverse Primer, Start at End of Read
				i = -(i+1)
				## Break After First 50 Bases
				if i <= -50:
					break
				## Mismatches and Ambiguities for Reverse
				for k in range(0,Rprimer_length):
					k = -(k+1)
					if RCr_primer[k] == "N":
						r_score = r_score + 1
					if RCr_primer[k] == "K":
						if read[i+(k+1)] == "A":
							r_score = r_score + 1	
						if read[i+(k+1)] == "C":
							r_score = r_score + 1
					if read[i+(k+1)] != RCr_primer[k]:
						r_score = r_score - 1
					## If Fall Below Threshold, Break
					if r_score < -(int(threshold)):
						break
					## Keep if Read Remains Above Threshold AND Spans Full Primer Length
					if k == -(Rprimer_length):
						rvs = 1
						## Cut Reverse Primer
						cut_read = cut_read[:((i-Rprimer_length)+1)]
						cut_read_length = len(cut_read)	
						## Check Length Requirement
						if cut_read_length < int(min_length):
							h = open("tempLog", "a+")
							h.writelines([key, "\nRead Did Not Pass Specified Length Requirement\n"])
							h.close()
							g = open(unk_file, "a+")
							g.writelines([key, "\n", read, "\n"])
							g.close()
						elif cut_read_length >= int(min_length):
							trim_count = trim_count + 1
							total_trimmed_length = total_trimmed_length + cut_read_length
							total_untrimmed_length = total_untrimmed_length + read_length
							f = open(output_file, "a+")
							f.writelines([key, "\n", cut_read, "\n"])
							f.close()
							break
		
		## If Forward Not Detected, Write to Log
		if fwd == 0:
			h = open("tempLog", "a+")
			h.writelines([key, "\nForward Primer Not Detected\n"])
			h.close()
			g = open(unk_file, "a+")
			g.writelines([key, "\n", read, "\n"])
			g.close()

		## If Forward Detected, but Not Reverse Write to Log
		elif rvs == 0:
			h = open("tempLog", "a+")
			h.writelines([key, "\nReverse Primer Not Detected\n"])
			h.close()
			g = open(unk_file, "a+")
			g.writelines([key, "\n", read, "\n"])
			g.close()

	
	## Statistics for LogFile
	e = open("mycut_jake-blamer_log.txt", "a+")
	
	## Length of Dictionary equals Number of Processed Reads
	e.writelines("Processed Reads: " + str(len(fa_dict)) + "\n")

	## trim_count variable Accounted for Trimmed Reads
	e.writelines("Trimmed Reads: " + str(trim_count) + "\n")

	## Adverage Untrimmed Read Length
	adv_untrimmed_length = total_untrimmed_length / trim_count
	e.writelines("Adverage Read Length -- Raw: " + str(adv_untrimmed_length) + "\n")

	## Adverage Trimmed Read Length
	adv_trimmed_length = total_trimmed_length / trim_count
	e.writelines("Adverage Read Length -- Trimmed: " + str(adv_trimmed_length) + "\n" + "Unprocessed Read Reasons: \n")

	## Combine Statistics with Explinations for Unprocessed Reads
	h = open("tempLog")
	e.write(h.read())
	h.close()
	e.close()

	## Remove temp Log File
	os.remove("tempLog")



## Help Message for Incorrect Commands
else:
	print("\n\nIncorrect Call! Try this Format:\n\n")
	print("python3 mycut_j-blamer.py –in_file <input.fasta> –out_file <output.fasta> –unk_file <unk.fasta> –n_mismatch <no. mistmatches> –min_len <minimum length> -forward <FORWARD_PRIMER> -reverse <REVERSE_PRIMER>\n")
	print("...or Enter 'python3 mycut_j-blamer.py --h' for Help\n")
