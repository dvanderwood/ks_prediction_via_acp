'''
After using HMMER on EMBL-EBI (https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch) with an ACP profile. 
Download the full length fasta file of the output and unzip it.

This script will use that fasta, a KS sequence simulated from an hmm profile, and the acp_hit_*_locations.txt output from top_acp_hits_finder.py to find the protein sequences
which contain the top ACP hits and to gather the KS sequences which follow these top hits.

usage: python ks_finder_v2.py fullseq_fasta simulated_ks acp_hit_*_locations.txt group_number_to_mark_seqs
'''

import sys, os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline

###Functions

#Function for blasting for a specified domain
def blastForDomain(domain, input_fasta, min_gene_length, e_value):
	blast_cline = NcbiblastpCommandline('blastp', query=domain, subject=input_fasta, outfmt=5, out= 'temp.xml')
	blast_cline()
	blast_record = next(NCBIXML.parse(open('temp.xml', 'r')))

	found_seq_dict = {}

	print('\n\tSignature for %s\n' % (domain.split('/')[-1]))
	E_VALUE_THRESH = e_value
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESH:
				if (hsp.sbjct_end - (hsp.sbjct_start-1)) > int(min_gene_length):
					print('\n**Alignment**')
					print('sequence:', alignment.title)
					print('e value:', hsp.expect)
					print(hsp.query[0:75] + '...')
					print(hsp.match[0:75] + '...')
					print(hsp.sbjct[0:75] + '...')
					print('Matched from position %i to %i in the input sequence' % (hsp.sbjct_start, hsp.sbjct_end))

					name = '>' + alignment.title.split(' ')[0] + ':' + str(hsp.sbjct_start) + '-' + str(hsp.sbjct_end)
					seq = hsp.sbjct

					found_seq_dict[name] = seq

	os.remove('temp.xml')
	
	


	return(found_seq_dict)

#Function to read a fasta containing multiple sequences.
def FastaToList(fasta_file):
	sequence = ''
	name = ''
	seqs = {}
	with open(fasta_file) as input_fasta:
		for line in input_fasta:
			if line.startswith('>'):
				if name == '':
					name = line.strip()[1:]
				else:
					seqs[name] = sequence
					name = line.strip()[1:]
					sequence = ''
			else:
				sequence += line.strip()
					
	return seqs

#Read a fasta that contains one sequence
def singleSequenceFastaToString(fasta_file): 
	input_seq = ''
	with open(fasta_file) as input_fasta:
		for line in input_fasta:
			if not line.startswith('>'):
				input_seq = input_seq + line.strip('\n')
	return input_seq

###Arguement clean-up
potential_hit_proteins = sys.argv[1]
simulated_ks = sys.argv[2]
acp_hit_locations = sys.argv[3]
group_number = sys.argv[4]

###Read Files
ks_query = singleSequenceFastaToString(simulated_ks)
#potential_hit_proteins_dict = FastaToList(potential_hit_proteins)

acp_hit_dict = {}
with open(acp_hit_locations, 'r') as acp_locations_file: #Gather only the ending location of the ACPs identified
	next(acp_locations_file)
	for line in acp_locations_file:
		line_split = line.strip().split('\t')
		locs = line_split[1:]
		loc_ends = set()
		for pos in locs:
			loc_ends.add(pos.split('-')[1])
		acp_hit_dict[line_split[0]] = loc_ends


###Blast for the KSs
found_seqs = blastForDomain(simulated_ks, potential_hit_proteins, 200, 1e-10)


###Check if the KS match up with a found ACP and store them
confirmed_ks = {}

for key in found_seqs:
	protein_name = key.split(':')[0][1:]
	location_start = key.split(':')[1].split('-')[0]
	if protein_name in acp_hit_dict:
		for loc in acp_hit_dict[protein_name]:
			if int(loc) + 100 > int(location_start) > int(loc):
				confirmed_ks[key] = found_seqs[key]

with open('ks_which_follow_acp_hits.fa', 'w') as output_file:
	for key in confirmed_ks:
		output_file.write(key + '|group_' + group_number + '\n')
		output_file.write(confirmed_ks[key] + '\n')



