#! /usr/bin/env python

import sys
import re
for line in sys.stdin:
    if not line.startswith('\t'):
        line_lst = line.strip().split('\t')
        temp = line_lst[0].split('|')
        if len(temp) == 10:  ## old version of AS output
            chrom, s1, s2, s3, s4, s5, s6, strand, name = temp[2],temp[4],temp[5],temp[6],temp[7],temp[8],temp[9],temp[3],temp[1]
        elif len(temp) == 7: ## new version of AS output
            chrom, s2, s3, s4, s5, strand, name = temp
            s1 = int(s2) - 100
            s6 = int(s5) + 100                
        s1, s2, s3, s4, s5, s6= [int(x) for x in [s1, s2, s3, s4, s5, s6]]
        change = line_lst[1]
        region1 = [s2-50, s2+300]
        region2 = [s3-300,s3+50]
        region3 = [s4-50,s4+300]
        region4 = [s5-300,s5+50]
        
        if s2-s1 < 50*2:
            region1[0] = s2-(s2-s1)/2
        if s3-s2 < 300*2:
            region1[1] = s2+(s3-s2)/2
            region2[0] = s3-(s3-s2)/2
        if s4-s3 < 50*2:
            region2[1] = s3+(s4-s3)/2
            region3[0] = s4-(s4-s3)/2
        if s5-s4 < 300*2:
            region3[1] = s4+(s5-s4)/2
            region4[0] = s5-(s5-s4)/2
        if s6-s5 < 50*2:
            region4[1] = s5+(s6-s5)/2
        
        step = 1
        
        for i, j in enumerate(range(s2,region1[0],-step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-step, j, line_lst[0], change, strand, 1, -i)
        for i, j in enumerate(range(s2,region1[1],step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+step, line_lst[0], change, strand, 1, i+1)

        for i, j in enumerate(range(s3,region2[0],-step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-step, j, line_lst[0], change, strand, 2, -i)
        for i, j in enumerate(range(s3,region2[1],step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+step, line_lst[0], change, strand, 2, i+1)
            
        for i, j in enumerate(range(s4,region3[0],-step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-step, j, line_lst[0], change, strand, 3, -i)
        for i, j in enumerate(range(s4,region3[1],step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+step, line_lst[0], change, strand, 3, i+1)
               
        for i, j in enumerate(range(s5,region4[0],-step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-step, j, line_lst[0], change, strand, 4, -i)
        for i, j in enumerate(range(s5,region4[1],step)):
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+step, line_lst[0], change, strand, 4, i+1)
            
