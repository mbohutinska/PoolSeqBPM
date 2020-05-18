################# C. AMARA SELECTION ANALYSIS BY MAJDA ######################
################# before start: prepare data ######################
setwd("/home/aa/Desktop/pracovni/cAmara/")
library(data.table)
pops<-c("CEZ","PIC","LUZ","VKR") #
for (p in pops) { # p="CEZ"
  a<-read.table(paste("AF/",p,"_AO.txt",sep=""),h=F)
  r<-read.table(paste("AF/",p,"_RO.txt",sep=""),h=F)
  i<-read.table(paste("AF/scaffold_position.txt",sep=""),h=F)
  all<-cbind(i,a,r)
  colnames(all)<-c('scaff','pos','a1','a2','a3','r1','r2','r3')
all<-subset(all,!all$a1 %like% ',')
all$af1<-as.numeric(as.character(all$a1))/(as.numeric(as.character(all$a1))+all$r1)  
all$af2<-as.numeric(as.character(all$a2))/(as.numeric(as.character(all$a2))+all$r2)  
all$af3<-as.numeric(as.character(all$a3))/(as.numeric(as.character(all$a3))+all$r3)  
t<-all[,9:11]
all$af<-apply(t,1,mean)
#this will change
#t$scaff<-paste("scaffold_1")
#t$pos<-1:nrow(t)
tot<-cbind(p,as.character(all$scaff),all$pos,all$af)
tot<-tot[complete.cases(tot), ]
write.table(x = tot,file = paste("AF/",p,"_AF.txt",sep=""),quote = F, row.names = F,col.names = F,sep = "\t")
}

################# 1. QUARTET FST ######################
### MAKE A SUMMARY MATRIX OF WINDOWS- 4 pops ###
setwd("/home/aa/Desktop/pracovni/cAmara/AF/")
library(data.table)
h1<-"CEZ"
h2<-"PIC"
f1<-"LUZ"
f2<-"VKR"
region<-"cAmara"

#read in all windows, make a list 
p1a<-fread(paste(h1,f1,"_WS1000_MS1_BPM.txt",sep=""),h=T)
p2a<-fread(paste(h1,f2,"_WS1000_MS1_BPM.txt",sep=""),h=T)
p3a<-fread(paste(f1,h2,"_WS1000_MS1_BPM.txt",sep=""),h=T) #changed to fit alphabet
p4a<-fread(paste(h2,f2,"_WS1000_MS1_BPM.txt",sep=""),h=T)
n1a<-fread(paste(h1,h2,"_WS1000_MS1_BPM.txt",sep=""),h=T)
n2a<-fread(paste(f1,f2,"_WS1000_MS1_BPM.txt",sep=""),h=T)

p1<-subset(p1a,p1a$num_snps>=20)
p2<-subset(p2a,p2a$num_snps>=20)
p3<-subset(p3a,p3a$num_snps>=20)
p4<-subset(p4a,p4a$num_snps>=20)
n1<-subset(n1a,n1a$num_snps>=20)
n2<-subset(n2a,n2a$num_snps>=20)

a<-rbind(p1[,2:4],p2[,2:4],p3[,2:4],p4[,2:4],n1[,2:4],n2[,2:4])
a1<-a[ order(a[,1], a[,2]), ]
a1<-subset(a1,!a1$scaff %in% "Genome")
a<-setDT(a1[!duplicated(a1[,c('scaff','start','end')]),])

# cds<-fread("lyrataCDSs.txt",h=T)
setkey(p1, scaff, start, end)
setkey(p2, scaff, start, end)
setkey(p3, scaff, start, end)
setkey(p4, scaff, start, end)
setkey(n1, scaff, start, end)
setkey(n2, scaff, start, end)

mp11<-foverlaps(a, p1, type = "within")
mp21<-foverlaps(a, p2, type = "within")
mp31<-foverlaps(a, p3, type = "within")
mp41<-foverlaps(a, p4, type = "within")
mn11<-foverlaps(a, n1, type = "within")
mn21<-foverlaps(a, n2, type = "within")

mp1<-mp11[!duplicated(mp11[,c('scaff','i.start','i.end')]),] 
mp2<-mp21[!duplicated(mp21[,c('scaff','i.start','i.end')]),] 
mp3<-mp31[!duplicated(mp31[,c('scaff','i.start','i.end')]),] 
mp4<-mp41[!duplicated(mp41[,c('scaff','i.start','i.end')]),] 
mn1<-mn11[!duplicated(mn11[,c('scaff','i.start','i.end')]),] 
mn2<-mn21[!duplicated(mn21[,c('scaff','i.start','i.end')]),] 

a$p1_Fst<-mp1$FstN
a$p2_Fst<-mp2$FstN
a$p3_Fst<-mp3$FstN
a$p4_Fst<-mp4$FstN
a$n1_Fst<-mn1$FstN
a$n2_Fst<-mn2$FstN

a$p1_snps<-mp1$num_snps
a$p2_snps<-mp2$num_snps
a$p3_snps<-mp3$num_snps
a$p4_snps<-mp4$num_snps
a$n1_snps<-mn1$num_snps
a$n2_snps<-mn2$num_snps

write.table(a,append = F,file = paste("QuartetSummary",region,".txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)

### RANK THE WINDOWS - 4 pops ###
setwd("/home/aa/Desktop/pracovni/cAmara/AF/")
library(data.table)
region<-"cAmara"
a<-fread(paste("QuartetSummary",region,".txt",sep=''))
mean(a$p1_snps,na.rm = T)
mean(a$p2_snps,na.rm = T)
mean(a$p3_snps,na.rm = T)
mean(a$p4_snps,na.rm = T)
mean(a$n1_snps,na.rm = T)
mean(a$n2_snps,na.rm = T)
a<-a[order(a$p1_Fst,decreasing = T),]
a$rank_p1<-seq(1,nrow(a),1)
a<-a[order(a$p2_Fst,decreasing = T),]
a$rank_p2<-seq(1,nrow(a),1)
a<-a[order(a$p3_Fst,decreasing = T),]
a$rank_p3<-seq(1,nrow(a),1)
a<-a[order(a$p4_Fst,decreasing = T),]
a$rank_p4<-seq(1,nrow(a),1)
a$sumRank<-a$rank_p1+a$rank_p2+a$rank_p3+a$rank_p4
a<-a[order(a$sumRank,decreasing = F),]
n1<-quantile(a$n1_Fst,0.99,na.rm = T)
n2<-quantile(a$n2_Fst,0.99,na.rm = T)
a$out_n1<-0
a$out_n1[a$n1_Fst > n1] <- 1
a$out_n2<-0
a$out_n2[a$n2_Fst > n2] <- 1
#Add annotation
cds<-fread("/home/aa/Desktop/pracovni/cAmara/gene-annotFile_MajdaModif.txt",h=T)
setkey(cds, scaff, start, end)
aaa<-foverlaps(a, cds, type="any")
aaa$start[is.na(aaa$start)]<-99
aaa$start [aaa$start  < aaa$i.start] <-aaa$i.start [aaa$start  < aaa$i.start]
aaa$end[is.na(aaa$end)]<--99
aaa$end [aaa$end  > aaa$i.end] <-aaa$i.end [aaa$end  > aaa$i.end]
aaa$sum<-aaa$end-(aaa$start)
aaa$sum[is.na(aaa$gene)]<-0
aa<-setDT(aaa)[,list(genic=sum(sum)),by=list(Category=paste(scaff,i.start))]
summary(aa$genic)
hist(aa$genic)
table (aa$genic)
a$genic<-aa$genic
#   b<-a[1:451,]
#   b1<-subset(b,b$out_n1 %in% 0 & b$out_n2 %in% 0)
a1<-subset(a,a$out_n1 %in% 0 & a$out_n2 %in% 0)
o1<-a1[1:(nrow(a1)*0.01),]
pdf(paste("genicVsFst_",region,".pdf",sep=""),width = 12,height = 8)
plot(a1$genic,a1$p1_Fst,ylab = "Fst",xlab="genoc space within 1000 bp window")
abline(h = min(o1$p1_Fst))
dev.off()
write.table(a,paste("Summary_rankFst",region,".txt",sep=""),quote=F,row.names = F,sep="\t")
#overlap with genes
o2<-foverlaps(o1, cds, type="any")
o3<-unique(o2$ID)
write.table(o3,paste("out_1percent_IDs_",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
write.table(o2,paste("out1percent_negOut_moreThan20SNPs_annotated",region,".txt",sep=''),quote=F,row.names = F,sep="\t")

### MAKE A SUMMARY MATRIX OF SNPS- 4 pops ###
setwd("/home/aa/Desktop/pracovni/cAmara/AF/")
library(data.table)
h1<-"CEZ"
h2<-"PIC"
f1<-"LUZ"
f2<-"VKR"
region<-"cAmara"

#read in all windows, make a list 
p1<-fread(paste(h1,f1,"_WS1_MS1_BPM.txt",sep=""),h=T)
p2<-fread(paste(h1,f2,"_WS1_MS1_BPM.txt",sep=""),h=T)
p3<-fread(paste(f1,h2,"_WS1_MS1_BPM.txt",sep=""),h=T)
p4<-fread(paste(h2,f2,"_WS1_MS1_BPM.txt",sep=""),h=T)
n1<-fread(paste(h1,h2,"_WS1_MS1_BPM.txt",sep=""),h=T)
n2<-fread(paste(f1,f2,"_WS1_MS1_BPM.txt",sep=""),h=T)

a<-rbind(p1[,2:4],p2[,2:4],p3[,2:4],p4[,2:4],n1[,2:4],n2[,2:4])
a1<-a[ order(a[,1], a[,2]), ]
a1<-subset(a1,!a1$scaff %in% "Genome")
a<-a1[!duplicated(a1[,c('scaff','start','end')]),]

# cds<-fread("lyrataCDSs.txt",h=T)
setkey(p1, scaff, start, end)
setkey(p2, scaff, start, end)
setkey(p3, scaff, start, end)
setkey(p4, scaff, start, end)
setkey(n1, scaff, start, end)
setkey(n2, scaff, start, end)

mp11<-foverlaps(a, p1, type = "within")
mp21<-foverlaps(a, p2, type = "within")
mp31<-foverlaps(a, p3, type = "within")
mp41<-foverlaps(a, p4, type = "within")
mn11<-foverlaps(a, n1, type = "within")
mn21<-foverlaps(a, n2, type = "within")

mp1<-mp11[!duplicated(mp11[,c('scaff','i.start','i.end')]),] 
mp2<-mp21[!duplicated(mp21[,c('scaff','i.start','i.end')]),] 
mp3<-mp31[!duplicated(mp31[,c('scaff','i.start','i.end')]),] 
mp4<-mp41[!duplicated(mp41[,c('scaff','i.start','i.end')]),] 
mn1<-mn11[!duplicated(mn11[,c('scaff','i.start','i.end')]),] 
mn2<-mn21[!duplicated(mn21[,c('scaff','i.start','i.end')]),] 

a$p1_Fst<-mp1$FstN
a$p2_Fst<-mp2$FstN
a$p3_Fst<-mp3$FstN
a$p4_Fst<-mp4$FstN
a$n1_Fst<-mn1$FstN
a$n2_Fst<-mn2$FstN

write.table(a,append = F,file = paste("QuartetSummarySNPs",region,".txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)

### RANK THE SNPs - 4 pops ###
setwd("/home/aa/Desktop/pracovni/cAmara/AF/")
library(data.table)
region<-"cAmara"
a<-fread(paste("QuartetSummarySNPs",region,".txt",sep=''))
a<-a[order(a$p1_Fst,decreasing = T),]
a$rank_p1<-seq(1,nrow(a),1)
a<-a[order(a$p2_Fst,decreasing = T),]
a$rank_p2<-seq(1,nrow(a),1)
a<-a[order(a$p3_Fst,decreasing = T),]
a$rank_p3<-seq(1,nrow(a),1)
a<-a[order(a$p4_Fst,decreasing = T),]
a$rank_p4<-seq(1,nrow(a),1)
a$sumRank<-a$rank_p1+a$rank_p2+a$rank_p3+a$rank_p4
a<-a[order(a$sumRank,decreasing = F),]
n1<-quantile(a$n1_Fst,0.99,na.rm = T)
n2<-quantile(a$n2_Fst,0.99,na.rm = T)
a$out_n1<-0
a$out_n1[a$n1_Fst > n1] <- 1
a$out_n2<-0
a$out_n2[a$n2_Fst > n2] <- 1
#Add annotation
cds<-fread("/home/aa/Desktop/pracovni/cAmara/gene-annotFile_MajdaModif.txt",h=T)
setkey(cds, scaff, start, end)
a1<-subset(a,a$out_n1 %in% 0 & a$out_n2 %in% 0)
o1<-a1[1:(nrow(a1)*0.01),]
write.table(a,paste("Summary_rankFstSNPs",region,".txt",sep=""),quote=F,row.names = F,sep="\t")

### OVERLAP WITH PIRITAS FINEMAV ###
fm<-fread("/home/aa/Desktop/pracovni/cAmara/annotationSummary_PiritaMAV_3.5.sorted.txt",h=F)
fm$scaff<-paste("scaffold_",fm$V1,sep="")
fm$start<-fm$V5
fm$end<-fm$V5
o1$start<-o1$end
setkey(fm, scaff, start, end)
o2<-foverlaps(o1, fm, type="any")
o3<-o2[complete.cases(o2[ , 2]),]
o3$one<-1
o4<-setDT(o3)[,list(perGene=sum(one)),by=list(Category=V6)]
o5<-subset(o4,o4$perGene>1)
o6<-subset(o4,o4$perGene>=3)

write.table(o3,paste("SNPs_Fst_MAV_Outliers_SortedByFst",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
write.table(o6,paste("topGenes3more_Fst_MAV_overlap",region,".txt",sep=''),quote=F,row.names = F,sep="\t")
write.table(o4,paste("topGenes1more_Fst_MAV_overlap",region,".txt",sep=''),quote=F,row.names = F,sep="\t")


### VennDiagram ###
library(VennDiagram)
fm$x<-paste(fm$scaff, fm$start,sep='.')
o1$x<-paste(o1$scaff, o1$start,sep='.')
#SNPs
v<-venn.diagram(x=list("SNP Fst"=o1$x,"MAV"=fm$x),paste("MAV_SNPFST_cAmara.tiff",sep=""),fill = c("gold","forestgreen"), lty = "blank",alpha=0.3 , cex=0.7, cat.cex=0.6,main = "Candidate SNPs",height = 2000,width = 2000)


###################### VISUALIZATION   ######################
# 1. candidate SNP-level - from fine-MAV - write Pirita?
## A. scatterplot: Fst
setwd("/home/aa/Desktop/pracovni/cAmara/plots/")
library(data.table)
library(ggplot2)
df<-fread("/home/aa/Desktop/pracovni/cAmara/AF/Summary_rankFstSNPscAmara.txt")

df<-df[ order(unlist(df[,1]),unlist(df[,2])), ]

an<-fread("/home/aa/Desktop/pracovni/cAmara/plots/orientationDec2019.txt")
fm<-fread("file:///home/aa/Desktop/pracovni/cAmara/annotationSummary_PiritaMAV_3.5.sorted.txt")

g<-read.table("/home/aa/Desktop/pracovni/cAmara/AF/out1percent_negOut_moreThan20SNPs_annotatedcAmara.txt",h=T)
genes<-g[!duplicated(g[,c(5)]),][,5]
genes<-na.omit(genes)

genes<-c("g8676","g28436","g18368","g23958","g23859","g4025","g5301","g29925","g10307") #top 9 arenosa
genes<-c("g1320","g16017","g20563","g3635") # other arenosa

for (ge in genes) { # ge<-"g15101"
  an1<-subset(an,substr(gene,4,10) %in% ge)
  df1<-subset(df,scaff %in% paste("scaffold_",an1$scaff,sep="") & start > (an1$start - 15000)  & end < (an1$end + 15000)) #zoomOut: 35000
  an2<-subset(an,scaff %in%  substr(df1$scaff,10,15) & an$start >= df1$start[1] & an$end <= df1$end[nrow(df1)])
  fm1<-subset(fm,fm$V4 %in% ge)
  df2<-df1[,4:7]
  df1$mean<-apply(df2,1,mean)
  dfm<-subset(df1,df1$end %in% fm1$V5)
  an2$mean<-max(as.numeric(as.character(na.omit(df1$mean))))
  if (nrow(an2)>0) {
    if (an1$orient %in% "+") {
      a<-ggplot() + 
        geom_point(data = df1, aes(x=as.numeric(as.character(end)), y=as.numeric(as.character(mean))),size = 3, alpha = 0.5) +
        geom_point(data = dfm, aes(x=as.numeric(as.character(end)), y=as.numeric(as.character(mean))),size = 3, alpha = 0.6,col="red3") + 
        geom_rect(aes(xmin = an1$start, xmax = an1$end, ymin = -Inf, ymax = Inf),alpha=0.1) +
        geom_segment(data = an2, aes(x = an2$start, y = an2$mean+0.02, xend = an2$end, yend = an2$mean+0.02), colour = "grey",size=1) +
        geom_segment(data = df1, aes(x = an1$start, y = max(as.numeric(as.character(mean)),na.rm = T)+0.02, xend = an1$end, yend = max(as.numeric(as.character(mean)),na.rm = T)+0.02), colour = "red3",size=1.5, arrow = arrow(length = unit(0.3, "cm"),ends = "last")) + 
        labs(x=paste("Position on scaffold ",an1$scaff,sep=""), y =expression('F'[ST]*' (diploids - tetraploids)'),title = paste(an1$ann)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size=12),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))
      
    } else
    { a<-ggplot() + 
      geom_point(data = df1, aes(x=as.numeric(as.character(end)), y=as.numeric(as.character(mean))),size = 3, alpha = 0.5) +
      geom_point(data = dfm, aes(x=as.numeric(as.character(end)), y=as.numeric(as.character(mean))),size = 3, alpha = 0.6,col="red3") +
      geom_rect(aes(xmin = an1$start, xmax = an1$end, ymin = -Inf, ymax = Inf),alpha=0.1) +
      geom_segment(data = an2, aes(x = an2$start, y = an2$mean+0.02, xend = an2$end, yend = an2$mean+0.02), colour = "grey",size=1) +
      geom_segment(data = df1, aes(x = an1$start, y = max(as.numeric(as.character(mean)),na.rm = T)+0.02, xend = an1$end, yend = max(as.numeric(as.character(mean)),na.rm = T)+0.02), colour = "red3",size=1.5, arrow = arrow(length = unit(0.3, "cm"),ends = "first")) + 
      labs(x=paste("Position on scaffold ",an1$scaff,sep=""), y =expression('F'[ST]*' (diploids - tetraploids)'),title = paste(an1$ann)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size=12),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))
    }
    
    pdf(paste(ge,".6Apr20.pdf",sep=""),height = 4,width = 8,pointsize = 12)
    print(a)
    dev.off()
  } else{}
}



## B. heatmap of AF
setwd("/home/aa/Desktop/pracovni/cAmara/plots/")
library(data.table)
library(stringr)
library(dplyr)
#Tetraploids
c<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/CEZ.table.recode.txt",h=F,na.strings = "-9")
c1<-c[,7:9]
c$mean<-apply(c1,1,mean)
c$af<-c$mean/c$V2
c1<-subset(c,c$mean > 20)
p<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/PIC.table.recode.txt",h=F,na.strings = "-9")
p1<-p[,7:9]
p$mean<-apply(p1,1,mean)
p$af<-p$mean/p$V2
p1<-subset(p,p$mean > 20)
v<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/VKR.table.recode.txt",h=F,na.strings = "-9")
v1<-v[,7:9]
v$mean<-apply(v1,1,mean)
v$af<-v$mean/v$V2
v1<-subset(v,v$mean > 20)
l<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/LUZ.table.recode.txt",h=F,na.strings = "-9")
l1<-l[,7:9]
l$mean<-apply(l1,1,mean)
l$af<-l$mean/l$V2
l1<-subset(l,l$mean > 20)
c2<-c[,c(3,4,12)]
p2<-p[,c(3,4,12)]
v2<-v[,c(3,4,12)]
l2<-l[,c(3,4,12)]
write.table(c2,"CEZaf.txt",row.names = F,col.names = c("scaff","start","af"))
write.table(p2,"PICaf.txt",row.names = F,col.names = c("scaff","start","af"))
write.table(v2,"VKRaf.txt",row.names = F,col.names = c("scaff","start","af"))
write.table(l2,"LUZaf.txt",row.names = F,col.names = c("scaff","start","af"))
### plot per gene - all sites in the gene
setwd("/home/aa/Desktop/pracovni/cAmara/plots/")
library(data.table)
library(gplots)
library(RColorBrewer)
c<-fread("CEZaf.txt",h=T)
p<-fread("PICaf.txt",h=T)
v<-fread("VRKaf.txt",h=T)
l<-fread("LUZaf.txt",h=T)
an<-fread("/home/aa/Desktop/pracovni/cAmara/gene-annotFile_MajdaModif.txt")

#SNPFst - fineMAV overlaps
genes<- c("g8401","g15174", "g9383","g8764","g1412","g4513","g9237", "g7530","g18050","g18056","g1400")
#All SNPFst
g<-read.table("genesSNPs.txt",h=F)
genes<-g[!duplicated(g[,c(1)]),]
#All WindowFst
g<-read.table("genesWindows.txt",h=F)
genes<-g[!duplicated(g[,c(1)]),]
#interesting WindowFst
genes<-c('g5732','g33987','g15174','g1430','g9237','g9237')

for (ge in genes) {
  an1<-subset(an,gene %in% ge)
  c1 <-subset(c,scaff %in% an1$scaff & start > an1$start  & start < an1$end)
  p1 <-subset(p,scaff %in% an1$scaff & start > an1$start  & start < an1$end)
  v1 <-subset(v,scaff %in% an1$scaff & start > an1$start  & start < an1$end)
  l1 <-subset(l,scaff %in% an1$scaff & start > an1$start  & start < an1$end)
  c1$end<-c1$start
  p1$end<-p1$start
  v1$end<-v1$start
  l1$end<-l1$start
  
  tt<-rbind(c1,p1,v1,l1)
  tt1<-tt[!duplicated(tt[,c(1,2)]),]
  setkey(c1, scaff, start, end)
  t0<-foverlaps(x = tt1, y = c1, type="within") #tt1, c1
  colnames(t0) <- c("scaff",seq(2,ncol(t0)-3,1),"start","af","end")
  setkey(p1, scaff, start, end)
  t1<-foverlaps(x = t0, y = p1, type="within") #c1, p1
  colnames(t1) <- c("scaff",seq(2,ncol(t1)-3,1),"start","af","end")
  setkey(v1, scaff, start, end)
  t2<-foverlaps(x = t1, y = v1, type="within") #c1, p1, v1
  colnames(t2) <- c("scaff",seq(2,ncol(t2)-3,1),"start","af","end")
  setkey(l1, scaff, start, end)
  t3<-foverlaps(x = t2, y = l1, type="within") #c1, p1, v1
  t3$tot<-(t3$`3`+t3$`9`+t3$`6`+t3$af)/4
    for (i in  1:nrow(t3)){ # i=1
      if (!t3$tot[i] %in% NA & t3$tot[i]>0.5)
      {t3[i,c(3,6,9,12)]<-1-t3[i,c(3,6,9,12)]
      } else {}}
  my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
  df<-as.matrix(t3[,c(12,9,6,3)],rownames = paste(t3$scaff,t3$start,sep=":"))
  pdf(paste("heatmap_",ge,"_",an1$ann,".pdf",sep=""),height = (nrow(t3)+10)/5,width = 7)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,2,4),sepcolor= c("black"),sepwidth = c (0.05),trace="none",ColSideColors = c("blue","blue","red","red"),labCol=c("CEZ","PIC","VKR","LUZ"),colCol= c("blue","blue","red","red"),offsetRow = c(0.05),offsetCol= c(1),cexRow = c(1),cexCol = c(2),lwid = c(0.05,0.8),margins = c(5,10),lhei = c(0.05,5),xlab = paste(an1$gene,an1$ann ,sep=" - "),srtCol=0)
  dev.off()
}

### plot per gene - only functional sites!!

setwd("/home/aa/Desktop/pracovni/cAmara/plots/")
library(data.table)
library(gplots)
library(RColorBrewer)
c<-fread("../AF/CEZ_AF.txt",h=F)
p<-fread("../AF/PIC_AF.txt",h=T)
v<-fread("../AF/VKR_AF.txt",h=T)
l<-fread("../AF/LUZ_AF.txt",h=T)
colnames(c)<-c("pop","scaff","start","af")
colnames(p)<-c("pop","scaff","start","af")
colnames(v)<-c("pop","scaff","start","af")
colnames(l)<-c("pop","scaff","start","af")
an<-fread("/home/aa/Desktop/pracovni/cAmara/plots/orientationDec2019.txt")

fun<-fread("/home/aa/Desktop/pracovni/cAmara/fineMAV/dip_tet_outliers_3.5_0.01.bin")
out<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/SNPs_Fst_MAV_Outliers_SortedByFstcAmara.txt")
fun1<-subset(fun, paste(fun$scaff,fun$pos,sep=".") %in% paste(out$V1,out$start,sep="."))
#1kb windows
g<-read.table("/home/aa/Desktop/pracovni/cAmara/AF/out1percent_negOut_moreThan20SNPs_annotatedcAmara.txt",h=T)
genes<-g[!duplicated(g[,c(5)]),][,5]

for (ge in genes) { # ge<-"g102"
  an1<-subset(an,substr(gene,4,10) %in% ge)
  c1 <-subset(c,scaff %in% paste("scaffold_",an1$scaff,sep="") & start > an1$start-2000  & start < an1$end+2000)
  p1 <-subset(p,scaff %in% paste("scaffold_",an1$scaff,sep="") & start > an1$start-2000  & start < an1$end+2000)
  v1 <-subset(v,scaff %in% paste("scaffold_",an1$scaff,sep="") & start > an1$start-2000  & start < an1$end+2000)
  l1 <-subset(l,scaff %in% paste("scaffold_",an1$scaff,sep="") & start > an1$start-2000  & start < an1$end+2000)
  
 # c1 <-subset(c1,paste(substr(c1$scaff,10,15),c1$start,sep=".") %in% paste(fun1$scaff,fun1$pos,sep="."))
 # p1 <-subset(p1,paste(substr(p1$scaff,10,15),p1$start,sep=".") %in% paste(fun1$scaff,fun1$pos,sep="."))
 # v1 <-subset(v1,paste(substr(v1$scaff,10,15),v1$start,sep=".") %in% paste(fun1$scaff,fun1$pos,sep="."))
 # l1 <-subset(l1,paste(substr(l1$scaff,10,15),l1$start,sep=".") %in% paste(fun1$scaff,fun1$pos,sep="."))
  c1$end<-c1$start
  p1$end<-p1$start
  v1$end<-v1$start
  l1$end<-l1$start
  
 # fun2 <-subset(fun1,paste(fun1$scaff,fun1$pos,sep=".") %in% paste(substr(c1$scaff,10,15),c1$start,sep="."))
 
 # if (nrow(fun2)>1) {
  tt<-rbind(c1,p1,v1,l1)
  tt1<-tt[!duplicated(tt[,c(1,2)]),]
  setkey(c1, scaff, start, end)
  t0<-foverlaps(x = tt1, y = c1, type="within") #tt1, c1
  colnames(t0) <- c("scaff",seq(2,ncol(t0)-3,1),"start","af","end")
  setkey(p1, scaff, start, end)
  t1<-foverlaps(x = t0, y = p1, type="within") #c1, p1
  colnames(t1) <- c("scaff",seq(2,ncol(t1)-3,1),"start","af","end")
  setkey(v1, scaff, start, end)
  t2<-foverlaps(x = t1, y = v1, type="within") #c1, p1, v1
  colnames(t2) <- c("scaff",seq(2,ncol(t2)-3,1),"start","af","end")
  setkey(l1, scaff, start, end)
  t3<-foverlaps(x = t2, y = l1, type="within") #c1, p1, v1
  #t3$AA<-substr(fun2$AA,4,15)
  t3$tot<-(t3$`3`+t3$`9`+t3$`6`+t3$af)/4
  t3$mav<-fun2$fineMAV
  t3$zero<-0
  for (i in  1:nrow(t3)){ # i=1
    if (!t3$tot[i] %in% NA & t3$tot[i]>0.5)
    {t3[i,c(3,6,9,12)]<-1-t3[i,c(3,6,9,12)]
    } else {}}
  my_palette <- colorRampPalette(c("khaki1", "green2", "blue3"))(n = 20)
  df<-as.matrix(t3[,c(12,9,6,3)],rownames = paste(t3$AA))
  pdf(paste("heatmapAAs_",ge,"_",an1$ann,".pdf",sep=""),height = (nrow(t3)+10)/5,width = 5)
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,2,4),sepcolor= c("black"),sepwidth = c (0.05),trace="none",ColSideColors = c("blue","blue","red","red"),labCol=c("CEZ","PIC","VKR","LUZ"),colCol= c("blue","blue","red","red"),offsetRow = c(0.05),offsetCol= c(1),cexRow = c(1),cexCol = c(2),lwid = c(0.05,0.8),margins = c(5,10),lhei = c(0.05,5),xlab = paste(an1$gene,an1$ann ,sep=" - "),srtCol=0)
  
 my_palette <- colorRampPalette(c("white","white","white","white","white", "yellow", "red"))(n = 40)
  df<-as.matrix(t3[,c(19,20)],rownames = round(t3$mav,digits = 1))
  heatmap.2(x = df,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,1),sepcolor= c("black"),sepwidth = c (0.05),trace="none",ColSideColors = c("white","white"),labCol=c("fineMAV",""),colCol= c("black"),offsetRow = c(0.05),offsetCol= c(1),cexRow = c(1),cexCol = c(2),lwid = c(0.05,0.8),margins = c(5,10),lhei = c(0.05,5),xlab = "",srtCol=0)
  dev.off()
  }  }


### 2. Manhattan Plot
setwd("/home/aa/Desktop/pracovni/cAmara/plots/")
library(data.table)
#tot<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/QuartetSummarycAmara.txt",h=T)
tot<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/QuartetSummarySNPscAmara.txt",h=T)
t1<-tot[,4:7]
t2<-tot[,8:9]
tot$FstP<-apply(t1,1,mean)  
tot$FstN<-apply(t2,1,mean)  
tot$Fst<-tot$FstP-tot$FstN
tot<-na.omit(tot)
tot<-subset(tot,tot$Fst>0.1)
hist(tot$Fst)
#install.packages("qqman")
library(qqman)
df<-tot[,c(3,1,2,12)]
df$scaff<-as.numeric(substr(df$scaff,10,15))
colnames(df)<-c("SNP","CHR","BP","P")
#df<-na.omit(df)
pdf("FstPerChrom1.pdf",width = 15,height = 3)
manhattan(df, chr="CHR", bp="BP", snp="SNP", p="P",logp = F,ylab="Fst")
dev.off()

df1<-subset(df,df$CHR < 200)
pdf("FstPerChrom2.pdf",width = 15,height = 3)
manhattan(df1, chr="CHR", bp="BP", snp="SNP", p="P",logp = F,ylab="Fst")
dev.off()

df1<-subset(df,df$CHR < 200)
df2<-subset(df1,df1$P > 0.25)

pdf("FstPerChrom3.pdf",width = 15,height = 3)
manhattan(df1, chr="CHR", bp="BP", snp="SNP", p="P",logp = F,ylab="delta Fst")
dev.off()

png("FstPerChrom3.png",width = 1920,height = 720,pointsize = 24)
manhattan(df2, chr="CHR", bp="BP", snp="SNP", p="P",logp = F,ylab="delta Fst",ylim=c(0.25,1.05))
dev.off()

df3<-subset(df2,df2$CHR < 91)
png("FstPerChrom4.png",width = 1920,height = 720,pointsize = 24)
manhattan(df3, chr="CHR", bp="BP", snp="SNP", p="P",logp = F,ylab="delta Fst",ylim=c(0.25,1.05))
dev.off()

table(df$CHR)

df3<-subset(df2,df2$CHR < 91)

cds<-fread("../gene-annotFile_MajdaModif.txt",h=T)
df4<-df3
df4$end<-df4$BP
df4$start<-df4$BP
df4$scaff<-paste("scaffold",df4$CHR,sep="_")
setkey(cds, scaff, start, end)
aaa<-foverlaps(df4, cds, type="any")
aaa$SNP<-aaa$gene
df4<-aaa[,c(7,8,9,10)]
png("FstPerChrom4.png",width = 1920,height = 720,pointsize = 24)
manhattan(df3, chr="CHR", bp="BP", snp="SNP", p="P",logp = F,ylab="delta Fst",ylim=c(0.25,1.05))
dev.off()



################# 1. QUARTET FST - INCL. DIVERSITY MEASURE ######################
### MAKE A SUMMARY MATRIX - 4 pops ###
setwd("/home/aa/alpine/scanTools/ScanTools/cAmara/")
library(data.table)
h1<-"CEZ"
h2<-"PIC"
f1<-"LUZ"
f2<-"VKR"
region<-"cAmara"

#read in all windows, make a list 
p1<-fread(paste(h1,f1,"_WS1000_MS10_BPM.txt",sep=""),h=T)
p2<-fread(paste(h1,f2,"_WS1000_MS10_BPM.txt",sep=""),h=T)
p3<-fread(paste(h2,f1,"_WS1000_MS10_BPM.txt",sep=""),h=T)
p4<-fread(paste(h2,f2,"_WS1000_MS10_BPM.txt",sep=""),h=T)
n1<-fread(paste(h1,h2,"_WS1000_MS10_BPM.txt",sep=""),h=T)
n2<-fread(paste(f1,f2,"_WS1000_MS10_BPM.txt",sep=""),h=T)
wf1<-fread(paste(f1,".WS1.0k_MS10_3ind_WPM.txt",sep=""),h=T)
wf2<-fread(paste(f2,".WS1.0k_MS10_3ind_WPM.txt",sep=""),h=T)
wh1<-fread(paste(h1,".WS1.0k_MS10_3ind_WPM.txt",sep=""),h=T)
wh2<-fread(paste(h2,".WS1.0k_MS10_3ind_WPM.txt",sep=""),h=T)
a<-rbind(p1[,2:4],p2[,2:4],p3[,2:4],p4[,2:4],n1[,2:4],n2[,2:4])
a1<-a[ order(a[,1], a[,2]), ]
a1<-subset(a1,!a1$scaff %in% "Genome")
a<-a1[!duplicated(a1[,c('scaff','start')]),] 

# cds<-fread("lyrataCDSs.txt",h=T)
setkey(p1, scaff, start, end)
setkey(p2, scaff, start, end)
setkey(p3, scaff, start, end)
setkey(p4, scaff, start, end)
setkey(n1, scaff, start, end)
setkey(n2, scaff, start, end)
setkey(wf1, scaff, start, end)
setkey(wf2, scaff, start, end)
setkey(wh1, scaff, start, end)
setkey(wh2, scaff, start, end)
setkey(a, scaff, start, end)
mp1<-foverlaps(a, p1, type="start")
mp2<-foverlaps(a, p2, type="start")
mp3<-foverlaps(a, p3, type="start")
mp4<-foverlaps(a, p4, type="start")
mn1<-foverlaps(a, n1, type="start")
mn2<-foverlaps(a, n2, type="start")
mwf1<-foverlaps(a, wf1, type="start",mult="first")
mwf2<-foverlaps(a, wf2, type="start",mult="first")
mwh1<-foverlaps(a, wh1, type="start",mult="first")
mwh2<-foverlaps(a, wh2, type="start",mult="first")
a$p1_Fst<-mp1$FstWC
a$p2_Fst<-mp2$FstWC
a$p3_Fst<-mp3$FstWC
a$p4_Fst<-mp4$FstWC
a$n1_Fst<-mn1$FstWC
a$n2_Fst<-mn2$FstWC
a$p1_mis<-mp1$num_sites
a$p2_mis<-mp2$num_sites
a$p3_mis<-mp3$num_sites
a$p4_mis<-mp4$num_sites
a$n1_mis<-mn1$num_sites
a$n2_mis<-mn2$num_sites
a$p1_snps<-mp1$num_snps
a$p2_snps<-mp2$num_snps
a$p3_snps<-mp3$num_snps
a$p4_snps<-mp4$num_snps
a$n1_snps<-mn1$num_snps
a$n2_snps<-mn2$num_snps
a$wf1Div<-mwf1$Diversity
a$wf2Div<-mwf2$Diversity
a$wh1Div<-mwh1$Diversity
a$wh2Div<-mwh2$Diversity
write.table(a,append = F,file = paste("QuartetSummary",region,".txt",sep=""),quote = F, sep = "\t",col.names = T,row.names = F)

pdf("snpDensity.fixed.pdf")
hist(a$p1_snps,breaks = 100)
hist(a$p2_snps,breaks = 100)
hist(a$p3_snps,breaks = 100)
hist(a$p4_snps,breaks = 100)
hist(a$n1_snps,breaks = 100)
hist(a$n2_snps,breaks = 100)
dev.off()

### CHECK FST DATA ###
region<-"cAmara"
a<-read.table(paste("QuartetSummary",region,".txt",sep=""),h=T)
a$sum1<-(as.numeric(a$wh1)+as.numeric(a$wf1))
a$sum2<-(as.numeric(a$wh1)+as.numeric(a$wf2))
a$sum3<-(as.numeric(a$wh2)+as.numeric(a$wf1))
a$sum4<-(as.numeric(a$wh2)+as.numeric(a$wf2))
pdf("checkFstDatacAmara.pdf",width = 12,height = 8)
plot(a$p1_Fst~a$p1_mis,ylab = "Fst - p1",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p2_Fst~a$p2_mis,ylab = "Fst - p2",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p3_Fst~a$p3_mis,ylab = "Fst - p3",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p4_Fst~a$p4_mis,ylab = "Fst - p4",xlab="number of sites/window",xlim=c(0,14000),ylim=c(-0.08,1))
plot(a$p1_Fst~a$p1_mis,ylab = "Fst - p1",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
plot(a$p2_Fst~a$p2_mis,ylab = "Fst - p2",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
plot(a$p3_Fst~a$p3_mis,ylab = "Fst - p3",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
plot(a$p4_Fst~a$p4_mis,ylab = "Fst - p4",xlab="number of sites/window - xlim=2000",xlim=c(0,2000),ylim=c(-0.08,1))
bpar<-par()
par(mfrow=c(2,1))
par(mar = c(4.1, 4.1, 0.1, 0.1))
plot(a$p1_Fst~a$wh1,ylab = "Fst - p1",xlab="nucleotide diversity in h1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p1_Fst~a$wf1,ylab = "Fst - p1",xlab="nucleotide diversity in f1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p2_Fst~a$wh1,ylab = "Fst - p2",xlab="nucleotide diversity in h1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p2_Fst~a$wf2,ylab = "Fst - p2",xlab="nucleotide diversity in f2",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p3_Fst~a$wh2,ylab = "Fst - p3",xlab="nucleotide diversity in h2",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p3_Fst~a$wf1,ylab = "Fst - p3",xlab="nucleotide diversity in f1",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p4_Fst~a$wh2,ylab = "Fst - p4",xlab="nucleotide diversity in h2",xlim=c(0,0.5),ylim=c(-0.08,1))
plot(a$p4_Fst~a$wf2,ylab = "Fst - p4",xlab="nucleotide diversity in f2",xlim=c(0,0.5),ylim=c(-0.08,1))
par(bpar)
m<- lm(a$p1_Fst~a$sum1)
plot(a$p1_Fst~a$sum1,ylab = "Fst - p1",xlab="sum od nucl. diversity h1+f1")
abline(m, col = "red")
m<- lm(a$p2_Fst~a$sum2)
plot(a$p2_Fst~a$sum2,ylab = "Fst - p2",xlab="sum od nucl. diversity h1+f2")
abline(m, col = "red")
m<- lm(a$p3_Fst~a$sum3)
plot(a$p3_Fst~a$sum3,ylab = "Fst - p3",xlab="sum od nucl. diversity h2+f1")
abline(m, col = "red")
m<- lm(a$p4_Fst~a$sum4)
plot(a$p4_Fst~a$sum4,ylab = "Fst - p4",xlab="sum od nucl. diversity h2+f2")
abline(m, col = "red")
dev.off()
pdf("histFstFagaras.pdf")
hist(a$p1_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$p2_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$p3_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$p4_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$n1_Fst,breaks = 100,xlim=c(-0.2,1))
hist(a$n2_Fst,breaks = 100,xlim=c(-0.2,1))
dev.off()

#SNP density
region<-"cAmara"
#a<-read.table(paste("Fst_SNPdensity",region,".txt",sep=""),h=T)
pdf(paste("Fst_SNPdensity",region,".pdf",sep=""),width = 12,height = 8)

plot(a$p1_Fst~a$p1_snps,ylab = "Fst - p1",xlab="number of SNPs per 1kb window")
plot(a$p2_Fst~a$p2_snps,ylab = "Fst - p2",xlab="number of SNPs per 1kb window")
plot(a$p3_Fst~a$p3_snps,ylab = "Fst - p3",xlab="number of SNPs per 1kb window")
plot(a$p4_Fst~a$p4_snps,ylab = "Fst - p4",xlab="number of SNPs per 1kb window")
plot(a$n1_Fst~a$n1_snps,ylab = "Fst - n1",xlab="number of SNPs per 1kb window")
plot(a$n2_Fst~a$n2_snps,ylab = "Fst - n2",xlab="number of SNPs per 1kb window")
dev.off()

#############OVERLAP WITH ARENOSA #######################
library(data.table)
#aa<-read.table("/home/aa/JICAutumn2016/finalAnalysis29Apr/results/posSelection/tetraploid/all_dip_tet_0.99_overlaps#GF.txt",h=T)
#ab<-as.data.frame(table(aa$ALcode))
#ab<-subset(ab,ab$Freq > 1)
#ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
#aa<-subset(ann, ann$AL %in% ab$Var1)
#aa<-subset(aa,!aa$AT %in% "nnn")
#aa<-aa[!duplicated(aa[,2]),] 
#write.table(aa,"arenosa_2more_FG.txt", row.names = F,quote = F)

aa<-fread("/home/aa/Desktop/pracovni/cAmara/candGenesATorthologs_arenosa_FstN_1kbwindows.txt",sep="\t",h=T)
ca<-read.table("/home/aa/Desktop/pracovni/cAmara/candidatesFinal.txt",h=T)

par1<-subset(ca,substr(at,1,9) %in% substr(aa$AT,1,9))
par2<-subset(aa, substr(aa$AT,1,9) %in% substr(ca$at,1,9))
par1<-par1[ order(par1[,2]), ]
par2<-par2[ order(par2[,2]), ]

par<-cbind(par1, par2)

##### GO ENRICHMENT TOPGO #####
#BiocManager::install("topGO")
#BiocManager::install("Rgraphviz")
setwd("/home/aa/Desktop/pracovni/cAmara/")
library("biomaRt")
library(topGO)
library(data.table)
#dict<-fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
#collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart",dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
GTOGO <- GTOGO[GTOGO$go_id != '',]
geneID2GO <- by(GTOGO$go_id,GTOGO$ensembl_gene_id,function(x) as.character(x))
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))

#amara
s<-read.table(paste("/home/aa/Desktop/pracovni/cAmara/candidatesFinal.txt",sep=""),h=T)
sel<-substr(s$at,1,9)

#arenosa
aa<-fread("/home/aa/Desktop/pracovni/cAmara/candGenesATorthologs_arenosa_FstN_1kbwindows.txt",sep="\t",h=T)
sel<-substr(aa$AT,1,9)


int.genes <- factor(as.integer(all.genes %in% sel))
names(int.genes) = all.genes

  #go.obj <- new("topGOdata", ontology='BP', allGenes = int.genes, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=5) ## BP, MF and CC, nodeSize - at least this number of genes
  go.obj <- new("topGOdata", ontology='BP', allGenes = int.genes, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=150) ## 
  
  #resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
  #resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
  resultsFe <- runTest(go.obj, algorithm = "elim", statistic = "fisher") #More conservative
  resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
  
  allRes <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100)
  allRes <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 100)
  
  
  a<-subset(allRes, allRes$Annotated>=5  & as.numeric(allRes$elimFisher) <=0.05) #& allRes$Annotated<=200
  #pdf(paste(lin,"/graph",lin,".pdf",sep=""),width = 14,height = 14)
  #showSigOfNodes(go.obj, score(resultsFe), firstSigNodes = 5, useInfo = 'all')
  write.table(file=paste("genes_genic_topgo_arenosa150.V2.txt",sep=""),a,sep="\t",row.names=F)
  #dev.off()

  
  ###VENN
  library(VennDiagram)
  v<-venn.diagram(x=list("C.amara"=substr(s$at,1,9),"A.arenosa"=substr(aa$AT,1,9)),paste("arenosaAmaraGenes.tiff",sep=""),fill = c("gold","forestgreen"), lty = "blank",alpha=0.3 , main = "Candidate genes")
  write.table(file=paste("genes_arenosa.txt",sep=""),sel,sep="\t",row.names=F)
  write.table(file=paste("genes_amara.txt",sep=""),sel,sep="\t",row.names=F)
  aaa<-subset(s,substr(at,1,9) %in% substr(aa$AT,1,9))
  
  a<-read.table("genes_genic_topgo_arenosa50_classic.txt",h=T)
  c<-read.table("genes_genic_topgo_amara50_classic.txt",h=T)
  v<-venn.diagram(x=list("C.amara"=c$GO.ID,"A.arenosa"=a$GO.ID),paste("arenosaAmaraGOs.tiff",sep=""),fill = c("gold","forestgreen"), lty = "blank",alpha=0.3 , main = "Candidate GOs")
  aa<-subset(a,GO.ID %in% c$GO.ID)
  
  a<-read.table("genes_genic_topgo_arenosa150.V2.txt",h=T)
  c<-read.table("genes_genic_topgo_amara150.txt",h=T)
 par<- subset(c,c$GO.ID%in%a$GO.ID)
  
  
  
  ####Signif. parallelism
 ##genes
 fisher.test(matrix(c(20000-6-229-452, 229-6, 452-6, 6), nrow=2), alternative="greater")
 
 ##GO
 fisher.test(matrix(c(6000-4-156-78, 78-4, 156-4, 4), nrow=2), alternative="greater")
 fisher.test(matrix(c(6000-5-73-22, 22-5, 73-5, 5), nrow=2), alternative="greater")
 
 
 
  p<- matrix(c(8, 794, 0, 10000), nrow = 2, dimnames = list(c("p", "no_p"), c("candidates", "no_candidates")))
  p
  fisher.test(p)
  
  p<- matrix(c(351, 2813-351, 2477, 2477517), nrow = 2, dimnames = list(c("p", "no_p"), c("candidates", "no_candidates")))
  p
  fisher.test(p)
  ###Overlap with fineMAV
  library(data.table)
  setwd("/home/aa/Desktop/pracovni/cAmara/")
  fm<-fread("/home/aa/Desktop/pracovni/cAmara/annotationSummary_PiritaMAV_3.5.sorted.txt",h=F)
  g<-read.table("amaraGeneList.txt")
  g$V2<-substr(g$V1,3,14)
  write.table("","candidates_fineMAVnubers.txt",col.names = F,row.names = F,quote = F)
  for (gg in g$V2) { ## gg="g11103"
    a<-subset(fm,fm$V4 %in% gg)
    write.table(paste(gg,nrow(a),sep="\t"),"candidates_fineMAVnubers.txt",col.names = F,row.names = F,quote = F,append = T)
    write.table(a,"candidates_fineMAVSNPs.txt",col.names = F,row.names = F,quote = F,append = T)
  }
  
  ###fineMAV numbers GW
  library(data.table)
  setwd("/home/aa/Desktop/pracovni/cAmara/")
  fm<-fread("/home/aa/Desktop/pracovni/cAmara/annotationSummary_PiritaMAV_3.5.sorted.txt",h=F)
  a<-as.data.frame(table(fm$V4))
aa<-subset(a,a$Freq > 2)  

  
    
#### to add numers to paper
  #data
library(data.table)
aa<-fread("file:///home/aa/Desktop/pracovni/cAmara/cAmara/CEZ.table.recode.txt")
mean(aa$V6)
nrow(aa)
aa<-fread("file:///home/aa/Desktop/pracovni/cAmara/cAmara/PIC.table.recode.txt")
mean(aa$V6)
nrow(aa)
aa<-fread("file:///home/aa/Desktop/pracovni/cAmara/cAmara/VKR.table.recode.txt")
mean(aa$V6)
nrow(aa)
aa<-fread("file:///home/aa/Desktop/pracovni/cAmara/cAmara/LUZ.table.recode.txt")
mean(aa$V6)
nrow(aa)

#bpm
setwd("/home/aa/Desktop/pracovni/cAmara/AF/")
a<-read.table("PICVKR_WS1000_MS1_BPM.txt",h=T)
fst<-mean(a$FstN,na.rm = T)
afd<-mean(abs(a$AFD),na.rm = T)
fix<-sum(a$FixedDiff,na.rm = T)
write.table(paste(afd,fix,fst),"bpm.txt",quote = F,col.names = F,row.names = F,append = T)
a<-read.table("LUZVKR_WS1000_MS1_BPM.txt",h=T)
fst<-mean(a$FstN,na.rm = T)
afd<-mean(abs(a$AFD),na.rm = T)
fix<-sum(a$FixedDiff,na.rm = T)
write.table(paste(afd,fix,fst),"bpm.txt",quote = F,col.names = F,row.names = F,append = T)
a<-read.table("LUZPIC_WS1000_MS1_BPM.txt",h=T)
fst<-mean(a$FstN,na.rm = T)
afd<-mean(abs(a$AFD),na.rm = T)
fix<-sum(a$FixedDiff,na.rm = T)
write.table(paste(afd,fix,fst),"bpm.txt",quote = F,col.names = F,row.names = F,append = T)
a<-read.table("CEZVKR_WS1000_MS1_BPM.txt",h=T)
fst<-mean(a$FstN,na.rm = T)
afd<-mean(abs(a$AFD),na.rm = T)
fix<-sum(a$FixedDiff,na.rm = T)
write.table(paste(afd,fix,fst),"bpm.txt",quote = F,col.names = F,row.names = F,append = T)
a<-read.table("CEZLUZ_WS1000_MS1_BPM.txt",h=T)
fst<-mean(a$FstN,na.rm = T)
afd<-mean(abs(a$AFD),na.rm = T)
fix<-sum(a$FixedDiff,na.rm = T)
write.table(paste(afd,fix,fst),"bpm.txt",quote = F,col.names = F,row.names = F,append = T)
a<-read.table("CEZPIC_WS1000_MS1_BPM.txt",h=T)
fst<-mean(a$FstN,na.rm = T)
afd<-mean(abs(a$AFD),na.rm = T)
fix<-sum(a$FixedDiff,na.rm = T)
write.table(paste(afd,fix,fst),"bpm.txt",quote = F,col.names = F,row.names = F,append = T)

###windows stats
a<-fread("file:///home/aa/Desktop/pracovni/cAmara/AF/Summary_rankFstcAmara.txt")

##find meiosis candidates from arenosa
library(data.table)
setwd("/home/aa/Desktop/pracovni/cAmara/")
cand<-c("AL1G10680","AL1G35730","AL2G25920","AL2G37810","AL4G46460","AL4G47570","AL6G15380","AL6G30890","AL8G25590","AL8G25600","AL8G26680")
cand<-c("AL1G26770", "AL1G62040", "AL5G35750", "AL6G10420")
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
aa<-subset(ann, ann$AL %in% cand)
sel<-substr(aa$AT,1,9)
### TODO find in orthogroup file, plot (AL-AT-CA code)
ort<-read.table("Orthogroups.txt",h=T,sep = "\t")

o<-""
for (i in 1:length(cand)) { #   i=1
  o<-rbind(o,subset(ort,ort$Alyrata_384_v2.1.protein %like% cand[i]))
}
fm<-fread("/home/aa/Desktop/pracovni/cAmara/annotationSummary_PiritaMAV_3.5.sorted.txt",h=F)
genes<-c("g8676","g28436","g18368","g23958","g23859","g4025","g5301","g29925","g10307")
mav<-subset(fm,fm$V4 %in% genes)


### Interactions STRING
setwd("/home/aa/Desktop/pracovni/cAmara/")
library(data.table)
#subset hits to CA AA
ca<-read.table("amaraCandList23Dec.txt")
a<-fread("candGenesATorthologs_arenosa_FstN_1kbwindows.txt",sep="\t",h=T)
aa<-as.data.frame(substr(a$AT,1,9))
colnames(aa)<-"V1"
write.table(aa,"arenosaCandList8Jan.txt",quote = F,col.names = F,row.names = F)
s<-read.table("string_interactions8Jan.tsv",h=T)
s$x<-substr(s$node1_external_id,6,14)
s$y<-substr(s$node2_external_id,6,14)
#for each line in CA, find protein with interactor in AA but not CA
for (i in 1:nrow(ca)) { #  i=3
  g<-as.character(droplevels(ca$V1[i]))
  s1<-subset(s,s$x==g | s$y==g)
  s2<-subset(s1,!s1$x %in% ca$V1[-i])
  s3<-subset(s2,!s2$y %in% ca$V1[-i] )
  print(g)
  print(nrow(s2))
  print(nrow(s3))
  if (nrow(s3)>0) {
    sx<-subset(s3,!x%in%g)
    sy<-subset(s3,!y%in%g)
    ss<-rbind(sx,sy)
    t<-paste(g,nrow(s3),paste(paste(sx$x, collapse = ','),paste(sy$y, collapse = ','),sep=","),sep="\t")
    write.table(t,"cAmaraIntArenosa_8Jan_V3.txt",row.names = F,col.names = F,quote = F,append = T)
  }
}
a<-read.table("cAmaraIntArenosa_8Jan_V3.txt",h=F)
a1<-subset(a,a$V2>1)
c<-read.table("amaraCandList23Dec_topATcodes.txt")
a2<-subset(a,a$V1 %in%c$V1)
a3<-subset(a2,a2$V2>1)
write.table(a1,"cAmaraIntArenosa_1more8Jan_V3.txt",quote = F,row.names = F,col.names = F)
write.table(a2,"cAmaraIntArenosa_topCand8Jan_V3.txt",quote = F,row.names = F,col.names = F)
write.table(a3,"cAmaraIntArenosa_topCand1more8Jan_V3.txt",quote = F,row.names = F,col.names = F)

####Permutations
#select random set of 1000 genes (feasible size for STRING)
library(dplyr)
setwd("/home/aa/Desktop/pracovni/cAmara/")
ort<-read.table("Orthogroups.txt",h=T,sep = "\t")
ort1<-sample_n(tbl = ort, size = 2000)
ort1<-ort1[,c(2,3,5)]
ort2<-na.exclude(ort1)
at<-as.data.frame(substr(ort2$TAIR10_pep_20110103_representative_gene_model,1,9))
write.table(at, "randomThanliana.txt",col.names = F,row.names = F,quote = F)
#search for assiciation among the "at" set of genes in STRING (miltiple proteins, thaliana), download .tsv, delete "#" at the beginning

## RANDOM DRAWS
library(dplyr)
s<-read.table("string_interactions_random.tsv",h=T)
tot<-read.table("randomThanliana.txt")
d<-""
df<-""
#start loop
for (n in 1:1000) { # n=1
ca<-sample_n(tbl = tot, size = 229,replace = F)
aa<-sample_n(tbl = tot, size = 452,replace = F)
s$x<-substr(s$node1_external_id,6,14)
s$y<-substr(s$node2_external_id,6,14)
#for each line in CA, find protein with interactor in AA but not CA
for (i in 1:nrow(ca)) { #  i=1
  g<-as.character(droplevels(ca$V1[i]))
  s1<-subset(s,s$x==g | s$y==g)
  s2<-subset(s1,!s1$x %in% ca$V1[-i] | !s1$y %in% ca$V1[-i] )
  if (nrow(s2)>0) {
    sx<-subset(s2,!x%in%g)
    sy<-subset(s2,!y%in%g)
    ss<-rbind(sx,sy)
    t<-paste(g,nrow(s2),paste(paste(sx$x, collapse = ','),paste(sy$y, collapse = ','),sep=","),sep="\t")
    write.table(t,"random.txt",row.names = F,col.names = F,quote = F,append = T)
  }
}
a<-read.table("random.txt",h=F)
a1<-subset(a,a$V2>1)
d<-c(d,nrow(a1))
df<-c(df,nrow(a))
file.remove("random.txt")
}
### p-value
d<-read.
distr<-as.numeric(df[2:length(df)])
observed = 104
  pdf("distributionAnyAssociation.pdf",width = 12,height = 4)
hist(distr,breaks = 20,xlim = c(0,observed))
abline(v=observed,col="red")
dev.off()
#calculate the pvalue (two sided test)
pvalue_observed_sampling <- function(distr, observed){
  pvalue <- 2 * min(sum(observed > distr) , sum(observed <= distr)) / length(distr)
  if(pvalue == 0 ) pvalue <- 1 / length(distr)
  return(pvalue)
}
pvalue_observed_sampling(distr, observed)
write.table(d,"distribution1More8Jan.txt")
write.table(df,"distribution0More8Jan.txt")


### Find diff SNPs for Sian
library(data.table)
#tot<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/QuartetSummarycAmara.txt",h=T)
tot<-fread("/home/aa/Desktop/pracovni/cAmara/cAmara/QuartetSummarySNPscAmara.txt",h=T)
ann<-fread("/home/aa/Desktop/pracovni/cAmara/dip_tet.table.txt",h=T)
or<-fread("file:///home/aa/Desktop/pracovni/cAmara/plots/orientationDec2019.txt")
cand<-fread("file:///home/aa/Desktop/pracovni/cAmara/topCandidates.txt",h=F)
or$ID<-substr(or$gene,4,15)
orcand<-subset(or,or$ID %in% cand$V1)
orcand$scaff<-paste("scaffold_",orcand$scaff,sep="")
ann<-ann[,1:7]
ann$scaff<-paste("scaffold_",ann$`#CHROM`,sep="")
ann$start<-ann$POS
ann$end<-ann$POS
t1<-tot[,4:7]
t2<-tot[,8:9]
tot$FstP<-apply(t1,1,mean)  
tot$FstN<-apply(t2,1,mean)  
tot$Fst<-tot$FstP-tot$FstN
tot<-na.omit(tot)
tot3<-subset(tot,tot$Fst>0.3)
tot3$start<-tot3$end
setkey(ann, scaff, start, end)
m<-foverlaps(tot3, ann, type="start")
maas<-subset(m,!m$POS %in% NA)
maas<-maas[,c(1:10,21)]
setkey(orcand, scaff, start, end)
mt<-foverlaps(maas, orcand, type="any")
mt1<-subset(mt,!mt$V2 %in% NA)
tot4<-subset(mt1,Fst>0.4)
tot5<-subset(mt1,Fst>0.5)

write.table(mt1,"/home/aa/Desktop/pracovni/cAmara/missenseSNPs_Fst_0.3.txt",quote = F,row.names = F)
write.table(tot4,"/home/aa/Desktop/pracovni/cAmara/missenseSNPs_Fst_0.4.txt",quote = F,row.names = F)
write.table(tot5,"/home/aa/Desktop/pracovni/cAmara/missenseSNPs_Fst_0.5.txt",quote = F,row.names = F)

#######areosa MDOC
a<-read.table("file:///home/aa/Desktop/pracovni/cAmara/Arenosa_mdoc.txt",h=T)
aa<-subset(a,a$Subsampleddataset %in% 1)
mean(aa$DOCmean)

#########Prepare fineMAV only candidates
library(data.table)
aa<-read.table("/home/aa/Desktop/pracovni/cAmara/annotationSummary_PiritaMAV_3.5.sorted.txt",sep = "\t")
ab<-as.data.frame(table(aa$V4))
ab<-subset(ab,ab$Freq > 1)

ort<-read.table("/home/aa/Desktop/pracovni/cAmara/Orthogroups.txt",h=T,sep = "\t")
o<-subset(ort, ort$Camara_protein %like% ab$Var1)

o<-""
for (i in 1:length(ab$Var1)) { #   i=1
  
  o<-rbind(o,subset(ort,ort$Camara_protein %like% ab$Var1[i]))
}


#########Scan arenosa

#    test.calcPairwisebpmAnn(recode_dir= "bordel/VCF_all300_DP8.M0.5", pops=['TET','DIP'], window_size=1000, min_snps=20, mem=1, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="3:59:00",scratch_gb=1, print1=False)
setwd("/home/aa/Desktop/pracovni/cAmara/")
library(data.table)
a<-read.table("/home/aa/Desktop/pracovni/cAmara/TETDIP_WS1000_MS20_BPM.txt",h=T,sep="\t")
a<-fread("file:///home/aa/JICAutumn2016/finalAnalysis29Apr/allBPM/TETDIP_WS1_MS1_BPM.txt",h=T,fill = T)

sum(a$FixedDiff)
mean(abs(na.omit(a$AFD)))

a<-a[ order(a[,21],decreasing = T), ]
a1<-setDT(a[1:(nrow(a)*0.01),])
g<-fread("file:///home/aa/Desktop/references/lyrataV2/genesIDsLyV2.gff",h=F)
colnames(g)<-c("scaff","V1","V2","start","end",colnames(g)[6:9])
setkey(a1, scaff, start, end)
tts<-foverlaps(g, a1, type = "any")
tts<-na.omit(tts)
#tts<-tts[!duplicated(tts[,c('V9')]),] 
tts<-tts[,10:29]
tts$ID<-substr(tts$V9,4,12)
s<-as.data.frame(table(tts$ID))
ann<-fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt",h=T,quote="")
aa<-subset(ann, ann$AL %in% s$Var1)
aa<-subset(aa,!aa$AT %in% "nnn")
aa<-aa[!duplicated(aa[,2]),] 
write.table(aa,"candGenesATorthologs_arenosa_FstN_1kbwindows.txt", row.names = F,quote = F,sep="\t")
write.table(tts,"candGenesWindowInfo_arenosa_FstN_1kbwindows.txt", row.names = F,quote = F,sep="\t")
