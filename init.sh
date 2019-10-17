source activate gmatic

if [ ! -d kmer ]; then
	mkdir kmer bg_seq doc
	touch kmer/SE_inclusion_all  kmer/SE_inclusion_dn  kmer/SE_inclusion_up
fi
