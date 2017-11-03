#! /usr/bin/env Rscript

library(RColorBrewer)
cols <- brewer.pal(3, 'Set1')
red <- cols[1]
blue <- cols[2]
green <- cols[3]
gray <- 'gray'
red2 <- rgb(col2rgb(cols[1])[1],col2rgb(cols[1])[2],col2rgb(cols[1])[3],150,max=255)
blue2 <- rgb(col2rgb(cols[2])[1],col2rgb(cols[2])[2],col2rgb(cols[2])[3],150,max=255)
#gray2 <- rgb(190,190,190,150,max=255)
gray2 <- 'black'

init <- function(file) {
  d <- read.table(file,row.names=1)
  x <- rep(0,350*4)
  names(x)  <- c(paste(1, seq(-49,300),sep='.'), paste(2, seq(-299,50),sep='.'), paste(3, seq(-49,300),sep='.'), paste(4, seq(-299,50),sep='.'))
  x[names(x) %in% rownames(d)] = d$V2
  return(x)
}

args <- commandArgs(T)
en_exon <- as.numeric(args[1])
si_exon <- as.numeric(args[2])
bg_exon <- as.numeric(args[3])
prefix <- '../RNAmap/'

# en_exon <- 345
# si_exon <- 186
# bg_exon <- 1963
# prefix <- './'

en_f <- init(paste0(prefix, "splice_sig_rnamap_en_f"))
en_r <- init(paste0(prefix, "splice_sig_rnamap_en_r"))
si_f <- init(paste0(prefix, "splice_sig_rnamap_si_f"))
si_r <- init(paste0(prefix, "splice_sig_rnamap_si_r"))
bg_f <- init(paste0(prefix, "splice_const_rnamap_f"))
bg_r <- init(paste0(prefix, "splice_const_rnamap_r"))
output <- paste0(prefix, "RNAmap_normByRNA.pdf")
output2 <- paste0(prefix, "RNAmap_3SS_normByRNA.pdf")

en_f_rna <- init(paste0(prefix, "splice_sig_rnamap_en_f.RNA"))
en_r_rna <- init(paste0(prefix, "splice_sig_rnamap_en_r.RNA"))
si_f_rna <- init(paste0(prefix, "splice_sig_rnamap_si_f.RNA"))
si_r_rna <- init(paste0(prefix, "splice_sig_rnamap_si_r.RNA"))
bg_f_rna <- init(paste0(prefix, "splice_const_rnamap_f.RNA"))
bg_r_rna <- init(paste0(prefix, "splice_const_rnamap_r.RNA"))

en <- (en_f + rev(en_r))
si <- (si_f + rev(si_r)) * -1
bg <- (bg_f + rev(bg_r))

en_rna <- (en_f_rna + rev(en_r_rna)) / 1000
si_rna <- (si_f_rna + rev(si_r_rna)) / 1000
bg_rna <- (bg_f_rna + rev(bg_r_rna)) / 1000

en <- en / en_rna
si <- si / si_rna
bg <- bg / bg_rna

n <- 350
gap <- 30
lwd <- 2

pdf(output,width=10,heigh=12)
par(mfrow=c(2,1))
par(mar=c(1,4,4,2))
y_min <- min(c(en, si))
y_max <- max(c(en, si))
y_range = c(y_min, y_max*1.4)

# alternative spliced exons
plot(-1000,-1000,xlab='',ylab='',ylim=y_range,xlim=c(0,(n+gap)*4),bty='n',xaxt='n',yaxt='n')
mtext('iclip tags/RNA-Seq reads*1000',side=2,at=y_min+diff(y_range)/2,line=2.5,cex=1.5)
# axis(2,at=pretty(c(0,y_max)),label=pretty(c(0,y_max)))
# axis(2,at=pretty(c(y_min,0)),label=-1*pretty(c(y_min,0)))
#axis(2,at=c(-80,-40,0,40,80,120),label=c(80,40,0,40,80,120),cex.axis=1.5)
axis(2,at=c(-3,-2,-1,0,1,2,3,4,5,6),label=c(3,2,1,0,1,2,3,4,5,6),cex.axis=1.5)
lines(1:n, en[1:n], col=red, lwd = lwd)
lines(1:n, si[1:n], col=blue, lwd = lwd)
lines((n+gap+1):(2*n+gap), en[(n+1):(n*2)], col=red, lwd = lwd)
lines((n+gap+1):(2*n+gap), si[(n+1):(n*2)], col=blue, lwd = lwd)
lines((2*n+2*gap+1):(3*n+2*gap), en[(n*2+1):(n*3)], col=red, lwd = lwd)
lines((2*n+2*gap+1):(3*n+2*gap), si[(n*2+1):(n*3)], col=blue, lwd = lwd)
lines((3*n+3*gap+1):(4*n+3*gap), en[(n*3+1):(n*4)], col=red, lwd = lwd)
lines((3*n+3*gap+1):(4*n+3*gap), si[(n*3+1):(n*4)], col=blue, lwd = lwd)

bars <- unique(c(seq(1,n,50), n, seq(n+gap,2*n+gap,50), 2*n+gap, seq(2*n+2*gap+1,3*n+2*gap,50), 3*n+2*gap, seq(3*n+3*gap+1,4*n+3*gap,50), 4*n+3*gap))
bars2 <- bars[c(1,2,4,6,8,9,11,13,15,16,17,18,20,22,24,25,27,29,31,32)] 
segments(bars2[1], 0, bars2[5], 0, lwd=lwd)
segments(bars2[6], 0, bars2[15], 0, lwd=lwd)
segments(bars2[16], 0, bars2[20], 0, lwd=lwd)
segments(bars2[5]-gap/2, y_max*0.05, bars2[5]+gap/2, -y_max*0.05, lwd=lwd)
segments(bars2[6]-gap/2, y_max*0.05, bars2[6]+gap/2, -y_max*0.05, lwd=lwd)
segments(bars2[15]-gap/2, y_max*0.05, bars2[15]+gap/2, -y_max*0.05, lwd=lwd)
segments(bars2[16]-gap/2, y_max*0.05, bars2[16]+gap/2, -y_max*0.05, lwd=lwd)
arrows(bars2[-c(5,6,15,16)], y_max*0.05, bars2[-c(5,6,15,16)], -y_max*0.05,length=0)
text(200,y_max,'100 nt',cex=1.5)
segments(150,y_max*0.9,250,y_max*0.9,lwd=2)
arrows(c(150,250),y_max*0.88,c(150,250),y_max*0.92,length=0)
s1 <- bars[1]
s2 <- bars[2]
s3 <- bars[15]
s4 <- bars[18]
s5 <- bars[31]
s6 <- bars[32]
segments(s2,y_max*1.17,s5,y_max*1.17, lwd=2)
polygon(c(s3,s3,s4),c(y_max*1.1,y_max*1.24,y_max*1.1),col=blue,border=NA)
polygon(c(s3,s4,s4),c(y_max*1.24,y_max*1.24,y_max*1.1),col=red,border=NA)
rect(s1,y_max*1.1,s2,y_max*1.24,col='gray',border=NA)
rect(s6,y_max*1.1,s5,y_max*1.24,col='gray',border=NA)
mtext(paste0('Enhanced alternative exons (n = ',en_exon,')'),side=3,col=red,line=0,cex=2)
mtext(paste0('Silenced alternative exons (n = ',si_exon,')'),side=1,col=blue,line=0,cex=2)

## Constitutively spliced exon
par(mar=c(3,4,0,2))
y_max <- max(bg)
y_range = c(0, y_max*1.4)
plot(-1000,-1000,xlab='',ylab='',ylim=y_range,xlim=c(0,(n+gap)*4),bty='n',xaxt='n',yaxt='n')
mtext('iclip tags/RNA-Seq reads*1000',side=2,at=y_max/2,line=2.5,cex=1.5)
axis(2,at=pretty(c(0,y_max)),label=pretty(c(0,y_max)),cex.axis=1.5)
#axis(2,at=pretty(c(y_min,0)),label=-1*pretty(c(y_min,0)))
lines(1:n, bg[1:n], lwd = lwd)
lines((n+gap+1):(2*n+gap), bg[(n+1):(n*2)], lwd = lwd)
lines((2*n+2*gap+1):(3*n+2*gap), bg[(n*2+1):(n*3)], lwd = lwd)
lines((3*n+3*gap+1):(4*n+3*gap), bg[(n*3+1):(n*4)], lwd = lwd)
bars2 <- bars[c(1,2,4,6,8,9,11,13,15,16,17,18,20,22,24,25,27,29,31,32)] 
segments(bars2[1], 0, bars2[5], 0, lwd=lwd)
segments(bars2[6], 0, bars2[15], 0, lwd=lwd)
segments(bars2[16], 0, bars2[20], 0, lwd=lwd)
segments(bars2[5]-gap/2, y_max*0.05, bars2[5]+gap/2, -y_max*0.05, lwd=lwd)
segments(bars2[6]-gap/2, y_max*0.05, bars2[6]+gap/2, -y_max*0.05, lwd=lwd)
segments(bars2[15]-gap/2, y_max*0.05, bars2[15]+gap/2, -y_max*0.05, lwd=lwd)
segments(bars2[16]-gap/2, y_max*0.05, bars2[16]+gap/2, -y_max*0.05, lwd=lwd)

arrows(bars2[-c(5,6,15,16)], y_max*0.05, bars2[-c(5,6,15,16)], 0,length=0)
text(200,y_max,'100 nt',cex=1.5)
segments(150,y_max*0.9,250,y_max*0.9,lwd=lwd)
arrows(c(150,250),y_max*0.88,c(150,250),y_max*0.92,length=0)
s1 <- bars[1]
s2 <- bars[2]
s3 <- bars[15]
s4 <- bars[18]
s5 <- bars[31]
s6 <- bars[32]
# segments(s2,y_min,s1,y_max, col=green)
# segments(s3,y_min,s2,y_max, col=green)
# segments(s4,y_min,s3,y_max, col=green)
# segments(s5,y_min,s4,y_max, col=green)
segments(s2,y_max*1.17,s5,y_max*1.17, lwd=2)
rect(s1,y_max*1.1,s2,y_max*1.24,col='black',border=NA)
rect(s3,y_max*1.1,s4,y_max*1.24,col='black',border=NA)
rect(s6,y_max*1.1,s5,y_max*1.24,col='black',border=NA)
segments((s3+s4)/2,y_max*1.1,(s3+s4)/2,y_max*1.24, lwd=2, col='white')
mtext(paste0('Constitutively spliced exon (n = ',bg_exon,')'),side=1,col=gray(0.2),line=1,cex=2)
dev.off()
# 
# pdf(output2,heigh=10,width=8)
# par(mfrow=c(2,1))
# plot(en[600:700],type='l',ylim=c(-20,20),col=red,lwd=lwd,ylab='crosslink site counts',yaxt='n',xaxt='n',xlab="3'SS",cex.lab=1.5)
# lines(si[600:700],col=blue,lwd=lwd)
# #plot(smooth.spline(600:700,en[600:700]),type='l',ylim=c(-20,20),col=red,lwd=lwd,ylab='crosslink site counts',yaxt='n',xaxt='n',xlab="3'SS",cex.lab=1.5)
# #lines(smooth.spline(600:700, si[600:700]),col=blue,lwd=lwd)
# abline(h=0)
# axis(2, c(-40,-20,0,20,40), c(40,20,0,20,40),cex.axis=1.5)
# axis(1, c(0,25,50,75,100), c(-50,-25,0,25,50),cex.axis=1.5)
# plot(bg[600:700],type='l',ylim=c(0,max(bg[600:700])),lwd=lwd,ylab='crosslink site counts',yaxt='n',xaxt='n',xlab="3'SS",cex.lab=1.5,cex.axis=1.5)
# #plot(smooth.spline(600:700,bg[600:700]),type='l',ylim=c(0,max(bg[600:700])),lwd=lwd,ylab='crosslink site counts',yaxt='n',xaxt='n',xlab="3'SS",cex.lab=1.5,cex.axis=1.5)
# abline(h=0)
# axis(1, c(0,25,50,75,100), c(-50,-25,0,25,50),cex.axis=1.5)
# axis(2, c(0,100,200,300), c(0,100,200,300),cex.axis=1.5)
# dev.off()
# 
