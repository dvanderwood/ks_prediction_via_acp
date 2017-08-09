'''
This script takes in the hmmsearch ouput file generated at the end of the ks_prediction_via_acp pipeline and generates the percent expected found KS versus KS found after ACPs from other groups above a cutoff bit score.

Currently the script uses a master list of suspcted mega groups and must be modified in the file.

Group numbers must be integers.

Usage: python3 hmmsearch_output_parser.py group_*_ks_seach.txt cutoff

group_*_ks_seach.txt = the output from hmmsearch
cutoff = an integer value of the bit score cutoff to compare number of expected found KS to unexpected

Limitations - Does not consider if the sequence was found under multiple UniProtKB entries or if the sequence was found by more than one ACP search. Both of these lead to the same sequence appearing
			  multiple times in the found_ks.afa file that was used as the seq. database for the hmmsearch.

'''

import sys, re

# Grab inputs and confirm the cutoff is an integer
hmmsearch_file = sys.argv[1]
if '/' in hmmsearch_file:
	file_iso = hmmsearch_file.split('/')[-1]
	group_number = int(file_iso.split('_')[1])
else:
	group_number = int(hmmsearch_file.split('_')[1])

cutoff = sys.argv[2]
try:
	int(cutoff)
except ValueError:
	sys.exit('\nError: The cutoff value is not an integer.\n')

cutoff = int(cutoff)

# Mega Groups List
mgroup = [
			[8,9,10],
			[2,11,25],
			[3,5,6]
		]
group_nums_to_check = set()
for g in mgroup:
	if group_number in g:
		group_nums_to_check = set(g)

if len(group_nums_to_check) < 1:
	group_nums_to_check = [group_number]


above_cutoff_hits = []
from_group_hits = 0
from_mgroup_hits = 0

# Read hmmsearch output file:
with open(hmmsearch_file,'r') as infile:
	for _ in range(15): #skip over unused initial lines in the output file, this assumes a consistent number of lines before the relevant information
		next(infile)

	for line in infile:
		if line == '\n': #stop reading the file after the score summary
			break


		line_split = line.strip().split() #split line by white characters
		
		if float(line_split[1]) > cutoff: #only look at hits over the cutoff for full sequence bit score
			seq_name = line_split[8]
			if len(seq_name.split('|')) <= 2: #ignore sequences used to construct any of the ACP profiles
				above_cutoff_hits.append(seq_name)
				
				if int(seq_name.split('|')[1].split('_')[1]) == group_number:
					from_group_hits += 1
				if len(group_nums_to_check) > 1:
					if int(seq_name.split('|')[1].split('_')[1]) in group_nums_to_check:
						from_mgroup_hits += 1
if len(above_cutoff_hits) > 0:
	percent_expected = from_group_hits/len(above_cutoff_hits) * 100
else:
	percent_expected = 0

print('\nOVERVIEW OF HITS FOR GROUP %s' % (group_number))
print('PERCENT HITS OVER %i FROM THE EXPECTED GROUP: %f' % (cutoff,percent_expected))
if len(group_nums_to_check) > 1: #if a the group is part of a mega group
	print('MEGA GROUP: ' + ', '.join(str(x) for x in group_nums_to_check))
	if len(above_cutoff_hits) > 0:
		percent_expected_mgroup = from_mgroup_hits/len(above_cutoff_hits) * 100
	else:
		percent_expected_mgroup = 0
	print('PERCENT HITS OVER %i FROM THE EXPECTED MEGA GROUP: %f' % (cutoff,percent_expected_mgroup))
print('TOTAL EXPECTED HITS OVER %i: %i' %(cutoff,from_group_hits))
print('TOTAL HITS OVER %i: %i\n' % (cutoff,len(above_cutoff_hits)))


