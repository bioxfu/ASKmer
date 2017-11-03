source activate gmatic
conda env export > doc/environment.yml

if [ ! -d kmer ]; then
	mkdir kmer bg_seq doc
fi
