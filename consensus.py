#!/usr/bin/env python

# Save as consensus.py, 
# run as python consensus.py input.fasta x
# where x is the percentage of sequences to call a position in the consensus sequence; 
# i.e. python consensus.py input.fasta 0.5 would mean that a residue or nucleotide would 
# have to be represented in 50% of the sequences to call that position.

import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

alignment = AlignIO.read(sys.argv[1], 'fasta')
summary_align = AlignInfo.SummaryInfo(alignment)
print summary_align.dumb_consensus(float(sys.argv[2]))
