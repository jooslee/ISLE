# Identifying clinically relevant Synthetic LEthality (ISLE)

# define relevant functions
# function for quantile normalization of gene expression and CNV
qnorm.array <- function(mat)
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}

# load libraries
library(data.table)
library(Rcpp)
library(survival)
library(parallel)

# enable openMP support for gcc
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

# set working directory
#setwd("/cbcb/project2-scratch/jooslee/SL/revision/github/")

# load TCGA data
load("data/prob.TCGA.RData")
genes=prob$genes
numGenes=length(prob$genes)
numSamples=length(prob$samples)

####### The mRNA expression (RNAseqV2) and patients' clinical characteristics were downloaded from TCGA Data Portal,
### and copy number data (Gistic values) was downloaded from Broad Institute's Firehose (https://gdac.broadinstitute.org/)
### on June 27, 2015. The TCGA data includes the following data fields:
### genes         19001 -none-     character: protein coding genes' symbols
### entrez        19001 -none-     numeric  : protein coding genes' entrez ID
### samples        8749 -none-     character: TCGA sample barcodes
### types          8749 -none-     character: cancer types 
### mRNA      166239749 -none-     numeric  : a matrix of gene expression (RNAseq) with genes in the row and samples in the column
### scna      166239749 -none-     numeric  : a matrix of SCNA (Gistic values) with genes in the row and samples in the column
### sex            8749 -none-     character: sex of patients
### age            8749 -none-     numeric  : age of patients 
### race           8749 -none-     character: ethnicity of patients
### survival      17498 -none-     numeric  : two column matrix with first column as survival time (in days) and second column as cencering (death=0, alive=1)
### stage          8749 -none-     character: stage of the tumor 
### mRNA.norm 166239749 -none-     numeric  : quantile-normalized of mRNA expression values
### scna.norm 166239749 -none-     numeric  : quantile-normalized of SCNA values
### surv.dt           2 data.table list		: survival matrix for Cox regression (first column: survival time (in days), second column: censoring (death=1, alive=0))
### GII            8749 -none-     numeric  : genomic instability index calculated following Bilal et al (2013).

### mRNAq2    166239749 -none-     numeric  : a matrix of mRNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### scnaq2    166239749 -none-     numeric  : a matrix of SCNA equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### for example, we used the following specific code to generate mRNAq2
###prob$mRNAq2=matrix(NA,dim(prob$mRNA)[1],dim(prob$mRNA)[2])
###utypes=unique(prob$types)
###for (i in seq(length(utypes))){
###	inx=which(prob$types==utypes[i])
###	mRNAtest=prob$mRNA[,inx]
###	mRNA1=mRNAtest> apply(mRNAtest, 1, quantile, probs = 1/3,  na.rm = TRUE)
###	mRNA2=mRNAtest> apply(mRNAtest, 1, quantile, probs = 2/3,  na.rm = TRUE)
###	prob$mRNAq2[,inx]=mRNA1+mRNA2
###	print(i)
###}

# define all possible SL interactions with all genes
sl=cbind(rep(seq(numGenes),each=numGenes),rep(seq(numGenes),numGenes))
ii=which(sl[,1]<sl[,2])
sl=sl[ii,]

# Step 0. define the initial pool
# (1) Initial Set I: Collect experimentally identified gold standard SLs
load("data/sl.golden.set.RData")
gd$sr=gd$sr[gd$flag==1,]
sl0=gd$sr
ix=which(sl0[,1]>sl0[,2])
sl0[ix,]=sl0[ix,c(2,1)]

# (2) Initial Set II: Infer candidate SLs from single gene knockout data
sourceCpp("src/essentialityTestPair.cpp",rebuild=T)
# col.names = (gene A SCNA, gene B ess) (gene A mRNA, gene B ess) 
# 			  (gene B SCNA, gene A ess) (gene B mRNA, gene A ess) for SL

load("data/prob.ach.old1.RData") # Cheung et al. PNAS (2011)
load("data/prob.ach.new1.RData") # Cowley et al. Sci. Data. (2014)
load("data/prob.ach.new2.RData") # Auirre et al. Cancer Discov (2016)
load("data/prob.mar.old1.RData") # Marcotte et al. Cancer Discov (2012)
load("data/prob.mar.new1.RData") # Marcotte et al. Cell (2016)
### We collected the following 5 large-scale shRNA/shRNA essentiality screening datasets for inferring the Initial Set II

### this table shows the data fields included in the first dataset (Cheung et al.) as an example,
### the remaining 4 datasets have a similar data structure.
###          Length  Class      Mode     
###genes       19001 -none-     character : gene symbols of protein coding genes
###entrez      19001 -none-     numeric   : entrez IDs of protein coding genes
###samples       102 -none-     character : sample IDs (as coded in the dataset)
###celllines     102 -none-     character : names of cell lines
###types         102 -none-     character : cancer types of cell lines
###mRNA      1938102 -none-     numeric   : matrix of mRNA expression with genes in rows, cell lines in columns
###mRNAq2    1938102 -none-     numeric   : matrix of mRNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across cell lines
###scna      1938102 -none-     numeric   : matrix of SCNA expression with genes in rows, cell lines in columns
###scnaq2    1938102 -none-     numeric   : matrix of SCNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across cell lines
###mut	     1938102 -none-     numeric   : matrix of mutation (1:mutation, NA:WT)
###mat           102 data.frame list      : matrix of shRNA/sgRNA essentiality scores with genes in rows, cell lines in columns (the lower the more essential)
FDR=0.2
pcl=NULL
for (k in seq(5)){
	if (k==1) probc=prob.ach.old1		
	if (k==2) probc=prob.ach.new1
	if (k==3) probc=prob.ach.new2
	if (k==4) probc=prob.mar.old1
	if (k==5) probc=prob.mar.new1
	mRNAq=probc$mRNAq2
	scnaq=probc$scnaq2
	mat=probc$mat
	p = essentialityTestPair(sl, t(scnaq), t(mRNAq), t(mat))
		# p[,1]: when gene 1 is more active, gene 2 has higher essentiality score (less essential) - SCNA
		# p[,2]: when gene 1 is more active, gene 2 has higher essentiality score (less essential) - mRNA
		# p[,3]: when gene 2 is more active, gene 1 has higher essentiality score (less essential) - SCNA
		# p[,4]: when gene 2 is more active, gene 1 has higher essentiality score (less essential) - mRNA
	pcl=cbind(pcl,p)
}
pcl[pcl< -999] = NA
pcl[pcl==0.5] = NA #when all values are the same, essentialityTestPair returns 0.5 p-value.

p.ess1=apply(pcl[,1:4],1,min,na.rm=T)
p.ess2=apply(pcl[,5:8],1,min,na.rm=T)
p.ess3=apply(pcl[,9:12],1,min,na.rm=T)
p.ess4=apply(pcl[,13:16],1,min,na.rm=T)
p.ess5=apply(pcl[,17:20],1,min,na.rm=T)
p.ess1[is.infinite(p.ess1)]=NA
p.ess2[is.infinite(p.ess2)]=NA
p.ess3[is.infinite(p.ess3)]=NA
p.ess4[is.infinite(p.ess4)]=NA
p.ess5[is.infinite(p.ess5)]=NA
	p.ess1.adj=p.adjust(p.ess1,"BH")
	p.ess2.adj=p.adjust(p.ess2,"BH")
	p.ess3.adj=p.adjust(p.ess3,"BH")
	p.ess4.adj=p.adjust(p.ess4,"BH")			
	p.ess5.adj=p.adjust(p.ess5,"BH")
	i.ess=which(p.ess1.adj<=FDR | p.ess2.adj<=FDR | p.ess3.adj<=FDR | p.ess4.adj<=FDR | p.ess5.adj<=FDR)
sl=sl[i.ess,]

# combine the initial pool from (1) and (2)
usr1=paste(sl[,1],sl[,2])
usr2=paste(sl0[,1],sl0[,2])
usr=union(usr1,usr2)
sl=do.call(rbind,strsplit(usr," "))
class(sl) <- "numeric"

# Step I. Identifying underrepresented pairs using hypergeometric test
sourceCpp("src/HyperGeometricTest.pair.cpp")
# molecular omics screening (based on mRNA)
mol.scna = hypergeometricTestPair(scnaq=prob$scnaq2, pairs=sl)
# molecular omics screening (based on SCNA)
mol.mRNA = hypergeometricTestPair(scnaq=prob$mRNAq2, pairs=sl)

# select candidate SL pairs that pass both mRNA or SCNA based test
i.mol=which(p.adjust(mol.scna[,1],"BH")<FDR & p.adjust(mol.mRNA[,1],"BH")<FDR)
sl=sl[i.mol,]

# Step II. Identifying survival-informative pairs using Cox regression
numGenes = nrow(prob$scna)
numSamples = ncol(prob$scna)
mRNA.norm=prob$mRNA.norm
scna.norm=prob$scna.norm
surv.dt=prob$surv.dt
age=qnorm.array(prob$age)
sex=as.character(prob$sex)
race=prob$race
types=prob$types
mRNAq2=prob$mRNAq2
scnaq2=prob$scnaq2
GII=qnorm.array(prob$GII)

cox.pair.sl = function(pair,use.mRNA=F)
{
        if(use.mRNA){
                g1 = mRNA.norm[pair[1],]
                g2 = mRNA.norm[pair[2],]
                f1 = mRNAq2[pair[1],]
                f2 = mRNAq2[pair[2],]
        }else{  
                g1 = scna.norm[pair[1],]
                g2 = scna.norm[pair[2],]
                f1 = scnaq2[pair[1],]
                f2 = scnaq2[pair[2],]
        }
        res=rep(c(0,1),2)
        if (sum(!is.na(g1))>100 & sum(!is.na(g2))>100 & sum(!is.na(f1))>100 & sum(!is.na(f2))>100){
        cov = ifelse(f1 ==0 & f2==0,1,0 )
        dt1 = data.frame(cbind(surv.dt, cov, age, sex, GII, race, types, cbind(g1 , g2)))
        cox.out = coxph(Surv(time,status) ~ cov + strata(types), data=dt1)
        aa  = summary(cox.out)        
        cox.out = coxph(Surv(time,status) ~ cov + g1 + g2 +strata(types), data=dt1)
        bb  = summary(cox.out)        
        cox.out = coxph(Surv(time,status) ~ cov + GII+age+ strata(types,sex,race), data=dt1)
        cc  = summary(cox.out)        
        cox.out = coxph(Surv(time,status) ~ cov + g1 + g2 +GII+age+strata(types,sex,race), data=dt1)
        dd  = summary(cox.out)        

        res=c(aa$coefficients["cov",c(1,5)], bb$coefficients["cov",c(1,5)],cc$coefficients["cov",c(1,5)], dd$coefficients["cov",c(1,5)])
	}
	return(res)
}

# survival screening (based on mRNA)
cox.mRNA = mclapply(1:nrow(sl), function(tt) cox.pair.sl(sl[tt,], use.mRNA=T),mc.cores=64)
cox.sl.mRNA = t(do.call(cbind, cox.mRNA))
# survival screening (based on SCNA)
cox.scna = mclapply(1:nrow(sl), function(tt) cox.pair.sl(sl[tt,], use.mRNA=F),mc.cores=64)
cox.sl.scna = t(do.call(cbind, cox.scna))

# select candidate SL pairs that pass either mRNA or SCNA based test
i.cox1=which(cox.sl.mRNA[,1]<0 & p.adjust(cox.sl.mRNA[,2],"BH")<FDR)
i.cox2=which(cox.sl.scna[,1]<0 & p.adjust(cox.sl.scna[,2],"BH")<FDR)
i.cox1=i.cox1[which(cox.sl.mRNA[i.cox1,7]<0 & p.adjust(cox.sl.mRNA[i.cox1,8],"BH")<FDR)]
i.cox2=i.cox2[which(cox.sl.scna[i.cox2,7]<0 & p.adjust(cox.sl.scna[i.cox2,8],"BH")<FDR)]
i.cox=union(i.cox1,i.cox2)
sl=sl[i.cox,]

# Step III: Identifying phylogeneticly linked pairs using NMF
load("data/yuval.phylogenetic.profile.RData")
### The phylogenetic profile is downloaded from Yuval Tabach et al. Mol Syst Biol. (2013), Supplementary Table 1
load("data/feature.weight.RData")
### the feature weights are determined based on the phylogenetic tree (Ensembl database: http://useast.ensembl.org/index.html)

### function to identify phylogenetic distance of a pair of genes
phylo.profile = function(sl.gene.all){
    sl.gene1 = sl.gene.all
    sl.phylo =  cbind(match(sl.gene1[,1], phylo$genes), match(sl.gene1[,2], phylo$genes))
    featureMat = (phylo[sl.phylo[,1],-(1:3)] - phylo[sl.phylo[,2],-(1:3)])^2
    featureMat %*% t(feature.weight)
}
pp=phylo.profile(cbind(prob$genes[sl[,1]],prob$genes[sl[,2]]))

# select candidate SL pairs that pass have small phylogenetic distance mesured by NMF (<10.5)
sl=sl[pp<10.5,] # this corresponds to about 50-percentile among the phylogenetic distance of all pairs

# Final clinical-SL network
cSL.pairs=cbind(prob$genes[sl[,1]],prob$genes[sl[,2]])


