#! /usr/bin/env python
#fromGTF.SE.bed
#ID  GeneID  geneSymbol  chr strand  exonStart_0base exonEnd upstreamES  upstreamEE  downstreamES    downstreamEE
#0   "ENSMUSG00000002015.5"  "Bcap31"    chrX    -   73712444    73712545    73706996    73707144    73714792    73714928
# R1: 5'-end of up-stream intron
# R2: 3'-end of up-stream intron
# R3: exon
# R4: 5'-end of down-stream intron
# R5: 3'-end of down-stream intron

import os
import sys

genome = sys.argv[1]
region_lens = [int(x) for x in sys.argv[4:]]
#region_len = 300
#region_len = 100

for region_len in region_lens:
    print(region_len)
    input_file = open(sys.argv[2])
    output_file = open(sys.argv[3]+'/SE_R1-R5.bed', 'w')

    for line in input_file:
        if line.startswith('ID'):
        	continue
        fields = line.strip('\n').split('\t')
        chrom, s1, s2, s3, s4, strand, gene = fields[3], fields[8], fields[5], fields[6], fields[9], fields[4], fields[2]
        gene = gene.replace('"', '')
        ID = chrom+"|"+s1+"|"+s2+"|"+s3+"|"+s4+"|"+strand+"|"+gene
        s1, s2, s3, s4 = [int(x) for x in [s1, s2, s3, s4]]
        if strand == '+':
            R1_start, R1_end = s1+10, s1+region_len
            R2_start, R2_end = s2-region_len, s2-31
            R3_start, R3_end = s2, s3 
            R4_start, R4_end = s3+10, s3+region_len
            R5_start, R5_end = s4-region_len, s4-31
        else:
            R5_start, R5_end = s1+10, s1+region_len
            R4_start, R4_end = s2-region_len, s2-31
            R3_start, R3_end = s2, s3
            R2_start, R2_end = s3+10, s3+region_len
            R1_start, R1_end = s4-region_len, s4-31
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, R1_start-1, R1_end, ID+':R1', '.', strand))
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, R2_start-1, R2_end, ID+':R2', '.', strand))
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, R3_start-1, R3_end, ID+':R3', '.', strand))
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, R4_start-1, R4_end, ID+':R4', '.', strand))
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, R5_start-1, R5_end, ID+':R5', '.', strand))
    output_file.close()
    input_file.close()

    os.system("bedtools getfasta -name -tab -s -fi "+genome+" -bed "+sys.argv[3]+"/SE_R1-R5.bed -fo "+sys.argv[3]+"/SE_R1-R5.temp")
    os.system("cat "+sys.argv[3]+"/SE_R1-R5.temp|sed 's/:/\t/'|sed -r 's/:.+\t/\t/' > "+sys.argv[3]+"/SE_R1-R5_"+str(region_len)+".seq")
    os.system("rm "+sys.argv[3]+"/SE_R1-R5.temp")
    os.system("rm "+sys.argv[3]+"/SE_R1-R5.bed")

