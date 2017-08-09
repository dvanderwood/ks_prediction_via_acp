'''
Use on the top_acp_*_hits.fa from top_acp_hits_finder.py along with the the unaligned fasta of 
ACPs used to build the hmm profile which was used to search UniProtKB to get the list of ACPs.

Gaps are removed during the fasta reading.

This will output a fasta identical to the first with a |original_name*** at the end of the names of the validated ACPs.

Usage: python acp_validator.py top_acp_*_hits.fa acps_which_built_profile
'''


import sys, os, shutil
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
#from Bio.Align.Applications import MuscleCommandline

###Functions

#Function for blasting for a specified domain
def blastForDomain_mark_exacts(domain, input_fasta):
	blast_cline = NcbiblastpCommandline('blastp', query=domain, subject=input_fasta, outfmt=5, out= 'temp.xml')
	blast_cline()
	blast_record = next(NCBIXML.parse(open('temp.xml', 'r')))

	validated = {}

	print('\n\tSignature for %s\n' % (domain.split('/')[-1]))
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.identities/hsp.align_length == 1:
				print('\n**Alignment**')
				print('sequence:', alignment.title)
				print('e-value:', hsp.expect)
				print(hsp.query[0:75] + '...')
				print(hsp.match[0:75] + '...')
				print(hsp.sbjct[0:75] + '...')
				print('identities:',hsp.identities)

				validated[alignment.title.split(' ')[0]] = domain.split('/')[-1].split('.')[0]

	os.remove('temp.xml')
	
	


	return(validated)

#Function to read a fasta containing multiple sequences and removes any gaps from the sequences
def FastaToList_gap_remover(fasta_file):
	sequence = ''
	name = ''
	seqs = {}
	with open(fasta_file) as input_fasta:
		for line in input_fasta:
			if line.startswith('>'):
				if name == '':
					name = line.strip()[1:]
				else:
					seqs[name] = sequence.replace('-','')
					name = line.strip()[1:]
					sequence = ''
			else:
				sequence += line.strip()

		seqs[name] = sequence.replace('-','')
					
	return seqs

###Arguement Clean-up
found_acps = sys.argv[1]
known_acps = sys.argv[2]

###Read Files
found_acps_dict = FastaToList_gap_remover(found_acps)
known_acps_dict = FastaToList_gap_remover(known_acps)

###Create temp files to use as blast queries
if not os.path.exists('temp_fasta'):
	os.makedirs('temp_fasta')

for query in known_acps_dict:
	with open('temp_fasta/' + query + '.fa', 'w') as temp_file:
		temp_file.write('>%s\n%s' % (query, known_acps_dict[query]))

###Run blast and search for only exact matches
exact_sequence_matches = {}
for key in known_acps_dict:
	 exact_sequence_matches_temp = blastForDomain_mark_exacts('temp_fasta/' + key + '.fa', found_acps)
	 exact_sequence_matches.update(exact_sequence_matches_temp)

shutil.rmtree('temp_fasta/')


###Show and check if all ACPs are found, if not show which were not found
acps_found = set()
#print('\nSequences with 1 to 1 identity to ACPs from the profile build stage\n')
for key in exact_sequence_matches:
	#print(key,exact_sequence_matches[key])
	acps_found.add(exact_sequence_matches[key])
if len(known_acps_dict) == len(acps_found):
	print('\nAll known ACPs found in %s' % found_acps)
else:
	to_del_list = []
	for known in known_acps_dict:
		if known in acps_found:
			to_del_list.append(known)
	print('\nThe following ACPs were not found in %s' % found_acps)
	for y in to_del_list:
		del known_acps_dict[y]
	for x in known_acps_dict:
		print(x)


###Write output fasta identical to the input but those which are identical to a known ACP are noted
with open('.'.join(found_acps.split('.')[0:-1]) + '_known_acps_marked.fa', 'w') as output_file:
	for name in found_acps_dict:
		if name in exact_sequence_matches:
			output_file.write('>' + name + '|' + exact_sequence_matches[name] + '***\n')
		else:
			output_file.write('>' + name + '\n')

		output_file.write(found_acps_dict[name] + '\n')






