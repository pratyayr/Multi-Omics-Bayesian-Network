
###### R code for choosing cohesive triplets, quadruples and performing 
###   the Bayesian Network Analysis starting from cleaned and normalized data #####


rm(list = ls())
#Set working directory as the same directory of the source file

### Load required packages ###
library("ggplot2")
library("Cairo")
library("bnlearn")
library("gplots")
library("psych")
library("lattice")
library("gridExtra")
library("reshape2")

### Run the source file containing the user-defined functions ###
source("functions.R")


### Load the saved R workspace containing all data ###
load("alldata.RObj") #It contains the following objects
#####################
# genomat.alleles.sdp.short: Genotype matrix
# genedata.collapsed: Gene (mRNA) expression matrix
# data.cp.collapsed: normalized and batch effect corrected miRNA expression matrix
# vsd.mirna: variance stabilized miRNA expression matrix
# mirna.names: list containing the miRNA names
# genenames: vector containing the gene names
# did2011: DID data
# lda2006: LDA data
# lorr2005: LORR data
#####################
load("allmodels.RObj") #This contains all exhaustive models for 4 node Bayesian Networks


###### Step 1: Identification of cohesive triplets #######

### Run the pairwise association between SDP-miRNA-mRNA ###
out.vst=geneeqtl.fisher(vsd.mirna,genomat.alleles.sdp.short) #miRNA-SDP
out.geneeqtl=geneeqtl.fisher(genedata.collapsed,genomat.alleles.sdp.short) #mRNA-SDP
out.mrna.mirna=geneeqtl.fisher(data.cp.collapsed,genedata.collapsed) #miRNA-mRNA

### P-value threshold ###
pval.thresh=1e-3
num.genes=27123
num.mir=881
num.snps=1416

genesnpsel=which(out.geneeqtl<pval.thresh,arr.ind=T)

candid=numeric() #Will store the cohesive triplets
for (i in 1:num.mir){
  if (floor(i/100)==i/100){print(i)}
  for (j in 1:nrow(genesnpsel)){
    pval.mirsnp=out.vst[i,genesnpsel[j,2]]
    pval.mirgene=out.mrna.mirna[i,genesnpsel[j,1]]
    pval.genesnp=out.geneeqtl[genesnpsel[j,1],genesnpsel[j,2]]
    if (pval.mirgene < pval.thresh & pval.mirsnp < pval.thresh){
      candid=rbind(candid,c(i,genesnpsel[j,2],genesnpsel[j,1]))
    }
  }
}

colnames(candid)=c("mir","sdp","gene")
print(dim(candid)) #25570 rows each storing a triplet. 





###### Step 2: Identification of quadruples for each phenotype and then 
###### Step 3: Bayesian Network Analysis (BNA) for each quadruple separately followed by
###### Step 4: Building larger Bayesian Networks. ###################

###### Now we will work pheotype by phenotype ######





#############################################
################## DID ######################
#############################################

is.element(row.names(did2011),colnames(vsd.mirna)) #Not all are present in our data
common=intersect(row.names(did2011),colnames(vsd.mirna))
counts=vsd.mirna[,common] #881 x 33
genos=genomat.alleles.sdp.short[,common]
genes=genedata.collapsed[,common]
did.pheno=did2011[common,] #length 33 vector
names(did.pheno)=common

didgenecor=didsnpcor=didmircor=numeric(nrow(candid)) 
#To store the correlations of the triplets with the phenotype
for (i in 1:nrow(candid)){
  didgenecor[i]=cor.test(did.pheno,genes[candid[i,3],])$p.value
  didsnpcor[i]=cor.test(did.pheno,as.numeric(genos[candid[i,2],]))$p.value
  didmircor[i]=cor.test(did.pheno,counts[candid[i,1],])$p.value
}

sum(didgenecor<0.05 & didsnpcor<0.05 & didmircor<0.05) #39
didcand=which(didgenecor<0.05 & didsnpcor<0.05 & didmircor<0.05)


###### Basic BNA for each quadruple separately ######
bnscoregap3=numeric(length(didcand))
for (i in 1:length(didcand)){
  temp=myhillclimb.exhaust(didcand[i],candid,counts,genos,genes,did.pheno,all.models,plot=F)
  bnscoregap3[i]=temp$bicgap[2]
}

sum(bnscoregap3>2) #0, no further investigation will be done







#############################################
################## LDA ######################
#############################################

common=intersect(row.names(lda2006),colnames(vsd.mirna))
counts=vsd.mirna[,common] #881 x 57
genos=genomat.alleles.sdp.short[,common]
genes=genedata.collapsed[,common]
lda.pheno=lda2006[common,] #length 57 vector
names(lda.pheno)=common

ldagenecor=ldasnpcor=ldamircor=numeric(nrow(candid))
#To store the correlations of the triplets with the phenotype
for (i in 1:nrow(candid)){
  ldagenecor[i]=cor.test(lda.pheno,genes[candid[i,3],])$p.value
  ldasnpcor[i]=cor.test(lda.pheno,as.numeric(genos[candid[i,2],]))$p.value
  ldamircor[i]=cor.test(lda.pheno,counts[candid[i,1],])$p.value
}

sum(ldagenecor<0.05 & ldasnpcor<0.05 & ldamircor<0.05) #2231
ldacand=which(ldagenecor<0.05 & ldasnpcor<0.05 & ldamircor<0.05)


###### Basic BNA for each quadruple separately ######
bnscoregap3=bnscoregap4=numeric(length(ldacand))
for (i in 1:length(ldacand)){
  temp=myhillclimb.exhaust(ldacand[i],candid,counts,genos,genes,lda.pheno,all.models,plot=F)
  bnscoregap3[i]=temp$bicgap[2]
  bnscoregap4[i]=temp$bicgap[3]
}

### Selecting interesting quadruples and building larger network ###
sel.mir=unique(candid[ldacand[which(bnscoregap3>2)],1])
sel.snps=unique(candid[ldacand[which(bnscoregap3>2)],2])
sel.genes=unique(candid[ldacand[which(bnscoregap3>2 & bnscoregap4>-1)],3])

smalldata=as.data.frame(cbind(t(genos[sel.snps,]),t(genes[sel.genes,]),
                              (counts[sel.mir,]),lda.pheno))
names(smalldata)=c(paste0("snp",1:length(sel.snps)),paste0("gene",1:length(sel.genes)),
                   paste0("mir",1:length(sel.mir)),"phe")

blackmat=blackmatgen(length(sel.snps),length(sel.genes),length(sel.mir))
bn.hc=hc(smalldata, blacklist=blackmat)
graphviz.plot(bn.hc)

# Find out which SDPs, genes and miRNAs are in the network #
row.names(genedata.collapsed)[sel.genes]
genenames[sel.genes]
mirna.names$names[sel.mir]
sel.snps


#Bootstrap#
set.seed(1)
boot=boot.strength(smalldata,R=500,algorithm="hc",
                   algorithm.args = list(blacklist=blackmat))
boot[(boot$strength > 0.5) & (boot$direction >= 0.5), ]
avg.boot=averaged.network(boot[boot$direction >= 0.5,],threshold = 0.5)
graphviz.plot(avg.boot)


####### Prediction #########

#### gene to phenotype ####
(summary(smalldata$gene2)[5]-summary(smalldata$gene2)[2])*(-1679.886)/
  (summary(smalldata$phe)[5]-summary(smalldata$phe)[2])

#### mirna to phenotype ####
(summary(smalldata$mir1)[5]-summary(smalldata$mir1)[2])*(0.0500073)*(-1679.886)/
  (summary(smalldata$phe)[5]-summary(smalldata$phe)[2])

### For Prediction figure ###
bn.fit(avg.boot,smalldata)
ldaarr.g=10708.563-1679.886*summary(smalldata$gene2)[c(1,2,5,6)]
ldaarr.m=10708.563-1679.886*(5.8883102+0.0500073*summary(smalldata$mir1)[c(1,2,5,6)])










#############################################
################## LORR ######################
#############################################


common=intersect(row.names(lorr2005),colnames(vsd.mirna))
lorr.pheno=lorr2005[common,] #Only males, 58 dimensional vector
names(lorr.pheno)=common
counts=vsd.mirna[,common] #881 x 58
genos=genomat.alleles.sdp.short[,common]
genes=genedata.collapsed[,common]


lorrgenecor=lorrsnpcor=lorrmircor=numeric(nrow(candid))
#To store the correlations of the triplets with the phenotype
for (i in 1:nrow(candid)){
  lorrgenecor[i]=cor.test(lorr.pheno,genes[candid[i,3],])$p.value
  lorrsnpcor[i]=cor.test(lorr.pheno,as.numeric(genos[candid[i,2],]))$p.value
  lorrmircor[i]=cor.test(lorr.pheno,counts[candid[i,1],])$p.value
}

sum(lorrgenecor<0.05 & lorrsnpcor<0.05 & lorrmircor<0.05) #646
lorrcand=which(lorrgenecor<0.05 & lorrsnpcor<0.05 & lorrmircor<0.05)



###### Basic BNA for each quadruple separately ######
bnscoregap3=bnscoregap4=numeric(length(lorrcand))
for (i in 1:length(lorrcand)){
  temp=myhillclimb.exhaust(lorrcand[i],candid,counts,genos,genes,lorr.pheno,all.models,plot=F)
  bnscoregap3[i]=temp$bicgap[2]
  bnscoregap4[i]=temp$bicgap[3]
}



### Selecting interesting quadruples and building larger network ###
sel.mir=unique(candid[lorrcand[which(bnscoregap3>2)],1])
sel.snps=unique(candid[lorrcand[which(bnscoregap3>2)],2])
sel.genes=unique(candid[lorrcand[which(bnscoregap3>2 & bnscoregap4>-1)],3])

## By checking the SDP locations, we figured out that there are two different SDP regions.
## Therefore we will fit two bigger networks

### Network 1 ###
sel.mir=sel.mir[-2]
sel.snps=sel.snps[-5]
sel.genes=sel.genes[-8]

smalldata=as.data.frame(cbind(t(genos[sel.snps,]),t(genes[sel.genes,]),
                              (counts[sel.mir,]),lorr.pheno))
names(smalldata)=c(paste0("snp",1:length(sel.snps)),paste0("gene",1:length(sel.genes)),
                   paste0("mir",1:length(sel.mir)),"phe")

blackmat=blackmatgen(length(sel.snps),length(sel.genes),length(sel.mir))
bn.hc=hc(smalldata, blacklist=blackmat)
graphviz.plot(bn.hc)


# Find out which SDPs, genes and miRNAs are in the network #
row.names(genedata.collapsed)[sel.genes]
genenames[sel.genes]
mirna.names$names[sel.mir]
sel.snps


#Bootstrap#
set.seed(1)
boot=boot.strength(smalldata,R=500,algorithm="hc",
                   algorithm.args = list(blacklist=blackmat))
boot[(boot$strength > 0.5) & (boot$direction >= 0.5), ]
avg.boot=averaged.network(boot[boot$direction >= 0.5,],threshold = 0.5)
graphviz.plot(avg.boot)




### Network 2 ###

sel.mir=unique(candid[lorrcand[which(bnscoregap3>2)],1])
sel.snps=unique(candid[lorrcand[which(bnscoregap3>2)],2])
sel.genes=unique(candid[lorrcand[which(bnscoregap3>2 & bnscoregap4>-1)],3])
sel.mir=sel.mir[2]
sel.snps=sel.snps[5]
sel.genes=sel.genes[8]

smalldata=as.data.frame(cbind(t(genos[sel.snps,]),(genes[sel.genes,]),
                              (counts[sel.mir,]),lorr.pheno))
names(smalldata)=c(paste0("snp",1:length(sel.snps)),paste0("gene",1:length(sel.genes)),
                   paste0("mir",1:length(sel.mir)),"phe")

blackmat=blackmatgen(length(sel.snps),length(sel.genes),length(sel.mir))
bn.hc=hc(smalldata, blacklist=blackmat)
graphviz.plot(bn.hc)


# Find out which SDPs, genes and miRNAs are in the network #
row.names(genedata.collapsed)[sel.genes]
genenames[sel.genes]
mirna.names$names[sel.mir]
sel.snps


#Bootstrap#
set.seed(1)
boot=boot.strength(smalldata,R=500,algorithm="hc",
                   algorithm.args = list(blacklist=blackmat))
boot[(boot$strength > 0.5) & (boot$direction >= 0.5), ]
avg.boot=averaged.network(boot[boot$direction >= 0.5,],threshold = 0.5)
graphviz.plot(avg.boot)


####### Prediction #########

#### gene to phenotype ####
(summary(smalldata$gene1)[5]-summary(smalldata$gene1)[2])*(-207.3)/
  (summary(smalldata$phe)[5]-summary(smalldata$phe)[2])

#### mirna to phenotype ####
(summary(smalldata$mir1)[5]-summary(smalldata$mir1)[2])*(-0.02816)*(-207.3)/
  (summary(smalldata$phe)[5]-summary(smalldata$phe)[2])

### For Prediction figure ###
bn.fit(avg.boot,smalldata)
lorrarr.g=1737.5735-207.2776*summary(smalldata$gene1)[c(1,2,5,6)]
lorrarr.m=1737.5735-207.2776*(8.13411157-0.02816277*summary(smalldata$mir1)[c(1,2,5,6)])

########## Analysis ends #######



##### Figures (not all figures were generated using R) #####

####### Prediction Figure #########
Cairo(file="ldalorrbar.pdf",typ="pdf",dpi=80,height=490,width=800)
op=par(oma=c(3.5,0,1,0))
par(mfrow=c(1,2))
plot(c(1,7.5),c(-600,600),type="n",
     xlab="Difference in Expression",ylab="Change in LDA (cm)",
     xaxt="n",main="A")
axis(side=1,at=c(2.5,6.5),labels = c("2-Quartile","4-Quartile"),tick = F)
lines(c(2,2),c(0,-ldaarr.m[2]+ldaarr.m[3]),col="red",lwd=20,lend="butt")
lines(c(2.7,2.7),c(0,-ldaarr.g[2]+ldaarr.g[3]),col="blue",lwd=20,lend="butt")
lines(c(6.1,6.1),c(0,-ldaarr.m[1]+ldaarr.m[4]),col="red",lwd=20, lend="butt")
lines(c(6.8,6.8),c(0,-ldaarr.g[1]+ldaarr.g[4]),col="blue",lwd=20, lend="butt")
abline(h=0)

plot(c(1,7.5),c(-50,50),type="n",
     xlab="Difference in Expression",ylab="Change in LORR (min)",
     xaxt="n",main="B")
axis(side=1,at=c(2.5,6.5),labels = c("2-Quartile","4-Quartile"),tick = F)
lines(c(2,2),c(0,-lorrarr.m[2]+lorrarr.m[3]),col="red",lwd=20, lend="butt")
lines(c(2.7,2.7),c(0,-lorrarr.g[2]+lorrarr.g[3]),col="blue",lwd=20, lend="butt")
lines(c(6.1,6.1),c(0,-lorrarr.m[1]+lorrarr.m[4]),col="red",lwd=20, lend="butt")
lines(c(6.8,6.8),c(0,-lorrarr.g[1]+lorrarr.g[4]),col="blue",lwd=20, lend="butt")
abline(h=0)
op = par(usr=c(0,1,0,1), xpd=NA)
legend(-0.5,-0.46,
       legend=c("miRNA","Gene"), 
       fill=c("red","blue"),
       bg="white", horiz=T, cex=1, bty="n")
dev.off()



#### Supplementary Figures 1-4 ####

#Correlation between the phenotypes
common=intersect(row.names(did2011),row.names(lda2006))
common1=intersect(common,row.names(lorr2005))

p1=did2011[common1,]
p2=lda2006[common1,]
p3=lorr2005[common1,]

Cairo(file="phenoscatter.pdf",typ="pdf",dpi=100,height=400,width=1200)
par(mfrow=c(1,3))
plot(p1,p2,pch=16,cex=0.8, col="darkgrey",xlab="DID",ylab="LDA",
     main=paste0("r = ",round(cor(p1,p2),2),
                 " (p = ",round(cor.test(p1,p2)$p.value,3),")"))
plot(p1,p3,pch=16,cex=0.8, col="darkgrey",xlab="DID",ylab="LORR",
     main=paste0("r = ",round(cor(p1,p3),2),
                 " (p = ",round(cor.test(p1,p3)$p.value,3),")"))
plot(p2,p3,pch=16,cex=0.8, col="darkgrey",xlab="LDA",ylab="LORR",
     main=paste0("r = ",round(cor(p3,p2),2),
                 " (p = ",round(cor.test(p3,p2)$p.value,3),")"))
dev.off()


### QQ plots ###
Cairo(file="QQplotconditional.pdf",typ="pdf",dpi=80,height=400,width=1200)
xx=out.vst[unique(which(out.mrna.mirna<0.001,arr.ind=T)[,1]),
           unique(which(out.geneeqtl<0.001,arr.ind=T)[,2])]
p1=qqunif.plot(as.numeric(xx),draw.conf=F,title="miRNA-SDP",yrange=c(0,10))
xx=out.geneeqtl[unique(which(out.mrna.mirna<0.001,arr.ind=T)[,2]),
                unique(which(out.vst<0.001,arr.ind=T)[,2])]
set.seed(123)                        
p2=qqunif.plot(as.numeric(xx),
               draw.conf=F,title="mRNA-SDP",yrange=c(0,10))
set.seed(123)
xx=out.mrna.mirna[unique(which(out.vst<0.001,arr.ind=T)[,1]),
                  unique(which(out.geneeqtl<0.001,arr.ind=T)[,1])]
p3=qqunif.plot(as.numeric(xx),
               draw.conf=F,title="miRNA-mRNA",yrange=c(0,10))
grid.arrange(p1,p2,p3,ncol=3)
dev.off()

Cairo(file="QQplotconditionalNZ.pdf",typ="pdf",dpi=80,height=400,width=1200)
xx=out.vst[unique(which(out.mrna.mirna<0.001,arr.ind=T)[,1]),
           unique(which(out.geneeqtl<0.001,arr.ind=T)[,2])]
p1=qqunif.plot(as.numeric(xx),draw.conf=F,title="miRNA-SDP")#,yrange=c(0,10))
xx=out.geneeqtl[unique(which(out.mrna.mirna<0.001,arr.ind=T)[,2]),
                unique(which(out.vst<0.001,arr.ind=T)[,2])]
set.seed(123)                        
p2=qqunif.plot(as.numeric(xx),
               draw.conf=F,title="mRNA-SDP")#,yrange=c(0,10))
set.seed(123)
xx=out.mrna.mirna[unique(which(out.vst<0.001,arr.ind=T)[,1]),
                  unique(which(out.geneeqtl<0.001,arr.ind=T)[,1])]
p3=qqunif.plot(as.numeric(xx),
               draw.conf=F,title="miRNA-mRNA")#,yrange=c(0,10))
grid.arrange(p1,p2,p3,ncol=3)
dev.off()


#Barplots#
dropstrains=c("ISS","ILS","D2","B6","DBA")
Cairo(file="strainbarplot.pdf",typ="pdf",dpi=100,height=1200,width=1400)
par(mfrow=c(3,1))
barplot(sort(did2011[-which(row.names(did2011) %in% dropstrains),1]),las=2,
        ylab="DID (g/kg)")
which(row.names(lorr2005) %in% dropstrains)#none
barplot(sort(lorr2005[,1]),las=2,ylab="LORR (min)")
which(row.names(lda2006) %in% dropstrains)#none
barplot(sort(lda2006[,1]),las=2,ylab="LDA (cm)")
dev.off()



### Histogram ###
Cairo(file="phenocorhist.pdf",typ="pdf",dpi=100,height=1200,width=1400)
par(mfrow=c(3,3))
#DID
common=intersect(row.names(did2011),colnames(vsd.mirna))
counts=vsd.mirna[,common] #881 x 33
genos=genomat.alleles.sdp.short[,common]
genes=genedata.collapsed[,common]
did.pheno=did2011[common,] #length 33 vector
names(did.pheno)=common


hist(cor(did.pheno,t(counts[unique(candid[,1]),])),main="DID-miRNA",xlab="Correlation")
hist(cor(did.pheno,t(genes[unique(candid[,3]),])),main="DID-mRNA",xlab="Correlation")
hist(cor(did.pheno,t(genos[unique(candid[,1]),])),main="DID-SDP",xlab="Correlation")


#LDA
common=intersect(row.names(lda2006),colnames(vsd.mirna))
counts=vsd.mirna[,common] #881 x 57
genos=genomat.alleles.sdp.short[,common]
genes=genedata.collapsed[,common]
lda.pheno=lda2006[common,] #length 57 vector
names(lda.pheno)=common

hist(cor(lda.pheno,t(counts[unique(candid[,1]),])),main="LDA-miRNA",xlab="Correlation")
hist(cor(lda.pheno,t(genes[unique(candid[,3]),])),main="LDA-mRNA",xlab="Correlation")
hist(cor(lda.pheno,t(genos[unique(candid[,2]),])),main="LDA-SDP",xlab="Correlation")


#LORR
common=intersect(row.names(lorr2005),colnames(vsd.mirna))
lorr.pheno=lorr2005[common,] #Only males, 58 dimensional vector
names(lorr.pheno)=common
counts=vsd.mirna[,common] #881 x 58
genos=genomat.alleles.sdp.short[,common]
genes=genedata.collapsed[,common]

hist(cor(lorr.pheno,t(counts[unique(candid[,1]),])),main="LORR-miRNA",xlab="Correlation")
hist(cor(lorr.pheno,t(genes[unique(candid[,3]),])),main="LORR-mRNA",xlab="Correlation")
hist(cor(lorr.pheno,t(genos[unique(candid[,2]),])),main="LORR-SDP",xlab="Correlation")

dev.off()


