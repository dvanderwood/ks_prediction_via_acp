# ks_prediction_via_acp
Trans-AT PKS KS domains can be predicted via ACP sequence.

Pipeline for Determining if KSs ACPs have a typical KS which follow them.

Required software: Python 3, Biopython, MUSCLE, HMMER3

Started with the ACPs and KSs for 26 modules (where the module is defined as processing domains, ACP, KS).

For each Group:
  - Use the default options hmmbuild (HMMER3) on both the alignments (MUSCLE) for the ACPs and KSs
  - Using hmmemit (HMMER3) and the -c option, create a consensus sequence using the KS HMM profile
  - Use https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch to search with the ACP HMM profile and download the XML and full length FASTA files
    - Options: UniProtKB for the database and Bit Score cutoffs
      - Cutoffs - Significance bit scores: Sequence = 60, Hit = 60
                - Report bit scores: Sequence = 40, Hit = 40
  - Run the top_acp_hits_finder.py on the XML file
    - A cutoff of 80 was used in this analysis
  - Run acp_validator.py on the FASTA output of top_acp_hits_finder.py and the ACPs which were used to construct the ACP HMM profile from before
    - In this case, those ACPs would be in the group_profile/group_#/acp_#.fa file
  - Use ks_finder_v2.py to find the KS domains on the same polypeptide as the found ACPs and that are within 100 amino acids after one of those same ACPs
    - Use the full length FASTA file, the consensus KS sequence, acp_hit_*_locations.txt output from the first script, and a group number identifier
  - Run profile_builder_check.py on the FASTA output of the last script and the KSs hich were used to construct the KS HMM profile from before
    - In this case, those KSs would be in the group_profile/group_#/ks_#.fa file

After the above is done for every group:
  - Concatanate all the output found and marked KS FASTAs into one and align (MUSCLE)
  - For each group, use hmmsearch (HMMER3) and the KS HMM profile so see the highest scoring found KSs for that profile


- Pull out KSs from sequences with ACP hits that are right after the ACP hits, use simulated KS from hmmemit - ks_finder_v2.py  - run profile_builder_check.py on the output to check for the KSs which would have been in our example modules
- Run KSs against ks hmm profiles made by hmmbuild using hmmsearch
- Check if over certain score

KS and ACP hmm profiles are built from the realigned fastas built for the preceding KS work

hmmbuild and hmmsearch all use standard parameters
hmmemit uses -c to build a consensus sequence


HMMER3 ONLINE DATA BASE GUIDELINES: 

UniProtKB

sig bit score: 60seq, 60 hit
report bit score: 40seq, 40 hit

So many hits the website was timing out when trying to download XML and Fullseq fasta with default settings

For KS search,

HMMER3 ONLINE DATA BASE GUIDELINES: for group_1 at least

UniProtKB

sig bit score: 400seq, 400 hit
report but score: 300seq, 300 hit
