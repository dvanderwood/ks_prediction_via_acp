'''
After using HMMER on EMBL-EBI (https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch) with an ACP profile. 
Download the xml file of the output. Given that xml file. This script will rank the individual hits via bit score.
Given a cutoff this script will then output the sequence containing the hit and the location of each hit within the sequence.


Usage: python top_acp_hits_finder.py -h
'''

import argparse, os

###Create arguement parser

parser = argparse.ArgumentParser(description='Find top scoring HMM sequence hits')

parser.add_argument('-p', '--path', type = str, help = 'Path where the xml file is located')

parser.add_argument('-c', '--cutoff', type = float, default = 60, help = 'The cutoff to decide to include in ranking. Default = 60.')

parser.add_argument('-o', '--output_location', type = str, default = os.getcwd(), help = 'Directory where the output file is written. Default is the current working directory.')

args = parser.parse_args()

###Arguement clean-up

xml_file = getattr(args, 'path')
score_cutoff = getattr(args, 'cutoff')
output_directory = getattr(args, 'output_location').rstrip('/')


class hmm_hit:
	'''Stores the the sequence hit information from the HMMER xml output'''

	def __init__(self, alisqname, iali, jali, aliaseq, bitscore, cevalue, ievalue):
		self.alisqname = alisqname
		self.iali = iali
		self.jali = jali
		self.aliaseq = aliaseq
		self.bitscore = bitscore
		self.cevalue = cevalue
		self.ievalue = ievalue

	def __repr__(self):
		return 'Protein: %s\nSequence Hit: %s\nStart: %s\nEnd: %s\nBitscore: %s\nce-value: %s\n' % (self.alisqname, self.aliaseq, self.iali, self.jali, self.bitscore, self.cevalue)



#proteins_containing_hits = {}
all_hits = []

with open(xml_file, 'r') as input_file:
	for line in input_file:
		if line.startswith('  <domains'): #Grab only lines showing hit data
			line_split = line.strip().split(' ')

			#Gather info need for hmm_hit class. This is very inneficient because different fields are contained for each xml hit. There is a better way to do this but that can be implemented later
			for field in line_split:
				if field.startswith('alisqname'):
					alisqname = field.split('=')[1].strip('"')
				elif field.startswith('iali'):
					iali = field.split('=')[1].strip('"')
				elif field.startswith('jali'):
					jali = field.split('=')[1].strip('"')
				elif field.startswith('aliaseq'):
					aliaseq = field.split('=')[1].strip('"')
				elif field.startswith('bitscore'):
					bitscore = field.split('=')[1].strip('"')
				elif field.startswith('cevalue'):
					cevalue = field.split('=')[1].strip('"')
				elif field.startswith('ievalue'):
					ievalue = field.split('=')[1].strip('"')

			#Debuging: Check if any field doesn' exist

			#Create hmm_hit class object for this hit and add it to the list

			hit = hmm_hit(alisqname,iali,jali,aliaseq,bitscore,cevalue,ievalue)

			#if not alisqname in proteins_containing_hits:
			#	proteins_containing_hits[alisqname] = [hit]
			#else:
			#	proteins_containing_hits[alisqname].append(hit)

			all_hits.append(hit)

#Sort all hits by bitscore
all_hits_sorted = sorted(all_hits, key = lambda x: float(x.bitscore), reverse = True)

#grab all hits above the cutoff and add their a location to a dictionary with the protein as the key
#Also ouput of all hits sorted by bitscore
hit_location_dict = {}

with open(output_directory + '/all_hits_sorted.txt', 'w') as all_hits_out_file:
	with open(output_directory + '/top_acp_' + str(score_cutoff) + '_hits.fa', 'w') as hits_fasta:
		for hit in all_hits_sorted:
			all_hits_out_file.write('%s: %s-%s\nSequence_Hit: %s\nBitscore: %s\n' % (hit.alisqname, hit.iali, hit.jali, hit.aliaseq, hit.bitscore))
			if float(hit.bitscore) > score_cutoff:
				hits_fasta.write('>%s|%s-%s\n%s\n' % (hit.alisqname, hit.iali, hit.jali, hit.aliaseq))
				if not hit.alisqname in hit_location_dict:
					hit_location_dict[hit.alisqname] = [hit.iali + '-' + hit.jali]
				else:
					hit_location_dict[hit.alisqname].append(hit.iali + '-' + hit.jali)


with open(output_directory + '/acp_hit_' + str(score_cutoff) + '_locations.txt', 'w') as output_file:
	output_file.write('Protein\tLocations\n')
	for protein in hit_location_dict:
		output_file.write(protein + '\t'+ '\t'.join(hit_location_dict[protein]) + '\n')

