# PoolSeqBPM
Scripts to calculate between-population genetic metrics from PoolSeq data

Author: M. Bohutinska

Scripts partly taken from P. Monnahans pipeline: https://github.com/pmonnahan/ScanTools

## 0. before start: prepare data in R
```
setwd("/home/aa/Desktop/pracovni/cAmara/")
library(data.table)
pops<-c("CEZ","PIC") #,"LUZ","VKR"
for (p in pops) { # p="CEZ"
  a<-read.table(paste("AF/",p,"_AO.txt",sep=""),h=F)
  r<-read.table(paste("AF/",p,"_RO.txt",sep=""),h=F)
  all<-cbind(a,r)
  colnames(all)<-c('a1','a2','a3','r1','r2','r3')
all<-subset(all,!all$a1 %like% ',')
all$af1<-as.numeric(as.character(all$a1))/(as.numeric(as.character(all$a1))+all$r1)  
all$af2<-as.numeric(as.character(all$a2))/(as.numeric(as.character(all$a2))+all$r2)  
all$af3<-as.numeric(as.character(all$a3))/(as.numeric(as.character(all$a3))+all$r3)  
t<-all[,7:9]
all$af<-apply(t,1,mean)
#this will change
t$scaff<-paste("scaffold_1")
t$pos<-1:nrow(t)
tot<-cbind(p,t$scaff,t$pos,all$af)
tot<-tot[complete.cases(tot), ]
write.table(x = tot,file = paste("AF/",p,"_AF.txt",sep=""),quote = F, row.names = F,col.names = F,sep = "\t")
}
```
You don't need to use this pre-processing script, just note that input for downstream scripts has 4 columns: pop name, scaffold, position and allele frequency (calculated as the average AF of all the pools; the AF in individual pools should be calculated as the fraction of total number of reads supporting the alternative allele).

## 1. combine per-population files with allele frequency together
```
sort -k2,2 -k3,3n -m CEZ_AF.txt PIC_AF.txt > CEZPIC.concat.txt
```
## 2. run the script
```
python3 bpm_PoolSeq.py -i CEZPIC.concat.txt -o ./ -prefix CEZPIC_WS1_MS1 -ws 1 -ms 1 -np 2
```
-i: input file generated in step 1.

-o: output directory

-prefix: name of the output file

-ws: window size in bp (i.e. 1000bp)

-ms: minimum number of SNPs required per window (i.e. 50)

-np: number of populations: 2

**Output: outname	scaff	start	end	win_size	num_snps	AFD	FixedDiff	FstN (Nei, 1987)


## 3. post-process
quartetScan.r summarises all scripts used to postprocess the Fst-based selection scans used in Bohutinska et al. 2020 https://doi.org/10.1101/2020.01.31.929109
