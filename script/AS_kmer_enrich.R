#! /usr/bin/env Rscript

library(Biostrings)

kmer_enrich = function(tab, id, region, k, permut_time) {
	
	seq_list = tab[tab$V1 %in% id & tab$V2 == region,]$V3
	kmer_mat = matrix(nrow = 4^k, ncol = length(seq_list))
	
	for (i in 1:length(seq_list)) {
		dna = seq_list[i]
		x = oligonucleotideFrequency(DNAString(dna),k)
		kmer_mat[,i] = x
		rownames(kmer_mat) = names(x)
	}
	kmer_sum = apply(kmer_mat, 1, sum)
	
	p = rep(0, 4^k)
	names(p) = names(kmer_sum)	
	
	step = permut_time/10
	for (n in 1:permut_time) {
		if (n %% (step) == 0) {
			cat(paste(k, ' mer ', region, ': ', n/step*10, "%\n", sep=''))
		}
		random_id = sample(tab$V1, length(id))
		random_list = tab[tab$V1 %in% random_id & tab$V2 == region,]$V3
		random_mat = matrix(nrow = 4^k, ncol = length(random_list))
		for (i in 1:length(random_list)) {
			random_dna=random_list[i]
			x = oligonucleotideFrequency(DNAString(random_dna),k)
			random_mat[,i] = x
			rownames(random_mat) = names(x)
		}
		random_sum = apply(random_mat, 1, sum)
		p[random_sum >= kmer_sum] = p[random_sum >= kmer_sum] + 1
	}
	list(kmer_sum, p/permut_time)
}

args = commandArgs(T)
AS_seq = args[1]
AS_lst = args[2]
output = args[3]

kmer = as.numeric(sub('k', '', strsplit(output, '.' ,fixed = T)[[1]][3]))
region_len = as.numeric(sub('bp', '', strsplit(output, '.' ,fixed = T)[[1]][2]))

permut = 1000

tab = read.table(AS_seq, stringsAsFactors=F)
ids = read.table(AS_lst)[,1]

Rs =  paste('R', 1:5, sep='')
dfm = data.frame(matrix(nrow=4^kmer, ncol=10, dimnames=list(names(oligonucleotideFrequency(DNAString('A'),kmer)), paste('R', rep(1:5,each=2), '_', rep(c('counts','pvalue'),5),sep=''))))

for (i in seq(1,10,by=2)) {
	x = kmer_enrich(tab, ids, Rs[(i+1)/2], kmer, permut)
	dfm[,i] = x[[1]]
	dfm[,i+1] = x[[2]]
}

write.table(dfm, output, quote=F, col.names=NA, sep='\t')


