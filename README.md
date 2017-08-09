# ks_prediction_via_acp
Trans-AT PKS KS domains can be predicted via ACP sequence.

Pipeline for Determining if KSs ACPs have a typical KS which follow them.

Required software: Python 3, Biopython, MUSCLE, HMMER3

Started with the ACPs and KSs for 26 modules (where the module is defined as processing domains, ACP, KS).

## For each Group:
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
  - Run profile_builder_check.py on the FASTA output of the last script and the KSs which were used to construct the KS HMM profile from before
    - In this case, those KSs would be in the group_profile/group_#/ks_#.fa file

## After the above is done for every group:
  - Concatenate all the output found and marked KS FASTAs into one and align (MUSCLE)
  - For each group, use hmmsearch (HMMER3) and the KS HMM profile so see the highest scoring found KSs for that profile
    - Each found KS is tagged by which group's ACP preceded it and the KSs with the same group number should be the highest scoring for that group's KS HMM profile
  - The output of each hmmsearch for each group was parsed with a hmmsearch_output_parser.py
    - KS hits with bit scores over 400 were selected
    - Any KS which were used to construct the HMM profile that performed that search were thrown out
    - The percent of these found KSs which came from the expected group was determined.
    - If a group was part of a suspected mega-group (See Final Notes), the percent of found KSs that were from any of the groups in that set was determined
    
## Final Notes
  - Some groups do not have a KS which follow the ACP or typically have KSs farther than 100 amino acids away or on another polypeptide and thus are not captured
    - In this analysis those groups are: 4, 11, 15, 18, 21, and 22
  - Some groups are very similar forming a set and thus the best matching KSs will come from any of the groups in the set
    - In this analysis those sets are
      - 8, 9, and 10
      - 2, 11, and 25
    - Individual members of group 7 are rather diverse and could be expected in the 2,11, and 25 set or possibly any of groups 3, 5, and 6
  - Some sequences are found by more than one ACP HMM profile search on Uniprot and are thus redundant but tagged by more than one group
    - This is seen by groups 3 and 5 for example
  - Some sequences are identical but have been submitted under different accession number and are not removed in this analysis
  - The last two points could be remedied in a more rigorous analysis
