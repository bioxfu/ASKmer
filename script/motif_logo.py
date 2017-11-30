import sys
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

input_file = sys.argv[1]
output_file = input_file+'.pdf'

aln = []
with open(input_file) as file:
	for line in file:
		aln.append(Seq(line.strip().replace('U','T')))

m = motifs.create(aln)
m.weblogo(output_file, format='pdf',  
	show_fineprint=False, show_ends=False, 
	show_errorbars=False, show_xaxis=False,
	yaxis_tic_interval=4)

