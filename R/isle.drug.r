# Identifying clinically relevant Synthetic LEthality (ISLE) for drug targets
# used for the patient drug response prediction
# query SL pairs: the [drug target] X [all genes] candidate SL pairs
# clinical screen was performed with the TCGA tumor samples w/o treatent record

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
#setwd("/cbcb/project2-scratch/jooslee/SL/revision/github3/")

# load TCGA data
load("data/prob.TCGA.RData")
genes=prob$genes
numGenes=length(prob$genes)
numSamples=length(prob$samples)
prob0=prob

####### The mRNA expression (RNAseqV2) and patients' clinical characteristics were downloaded from TCGA Data Portal,
### and copy number data (Gistic values) was downloaded from Broad Institute's Firehose (https://gdac.broadinstitute.org/)
### on June 27, 2015. The TCGA data includes the following data fields:
### genes         19001 -none-     character: protein coding genes' symbols
### entrez        19001 -none-     numeric  : protein coding genes' entrez ID
### samples        8749 -none-     character: TCGA sample barcodes
### types          8749 -none-     character: cancer types 
### mRNA      166239749 -none-     numeric  : a matrix of gene expression (RNAseq) with genes in the row and samples in the column
### scna      166239749 -none-     numeric  : a matrix of SCNA (Gistic values) with genes in the row and samples in the column
### mRNAq2    166239749 -none-     numeric  : a matrix of mRNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### scnaq2    166239749 -none-     numeric  : a matrix of SCNA equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### sex            8749 -none-     character: sex of patients
### age            8749 -none-     numeric  : age of patients 
### race           8749 -none-     character: ethnicity of patients
### survival      17498 -none-     numeric  : two column matrix with first column as survival time (in days) and second column as cencering (death=0, alive=1)
### stage          8749 -none-     character: stage of the tumor 
### mRNA.norm 166239749 -none-     numeric  : quantile-normalized of mRNA expression values
### scna.norm 166239749 -none-     numeric  : quantile-normalized of SCNA values
### surv.dt           2 data.table list		: survival matrix for Cox regression (first column: survival time (in days), second column: censoring (death=1, alive=0))
### GII            8749 -none-     numeric  : genomic instability index calculated following Bilal et al (2013).

load("data/TCGA.drug.response/prob.TCGA.drug.response.RData")
targetIDs=unique(prob$drug.target$targetIDs)
sl=cbind(rep(targetIDs,each=19001),rep(seq(19001),length(targetIDs)))
ix=which(sl[,1]!=sl[,2])
sl=sl[ix,]
	
# Step 0. define the initial pool
# (1) Initial Set I: Collect experimentally identified gold standard SL
prob=prob0
load("data/sl.golden.set.RData")
gd$sr=gd$sl[gd$flag==1,]
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
pcl[pcl< -999] = NA;pcl[pcl==0.5] = NA
p.ess=apply(cbind(pcl),1,min,na.rm=T)
	p.ess[is.infinite(p.ess)]=NA
	p.ess.adj=p.adjust(p.ess,"BH")
	
# combine the initial pool from (1) and (2)	
use1=c(paste(sl0[,1],sl0[,2]),paste(sl0[,2],sl0[,1]))
use2=paste(sl[,1],sl[,2])
iix=match(use1,use2)
iix=iix[!is.na(iix)]
	p.ess[iix]=1E-6
	p.ess.adj[iix]=1E-6
	i.ess=which(p.ess.adj<0.1)
sl.tot=cbind(sl,p.ess,p.ess.adj)
sl=sl.tot[,1:2]

# Step I. Identifying underrepresented pairs using hypergeometric test
sourceCpp("src/HyperGeometricTest.pair.cpp")
mol.scna = hypergeometricTestPair(scnaq=prob$scnaq2, pairs=sl)
mol.mRNA = hypergeometricTestPair(scnaq=prob$mRNAq2, pairs=sl)
sl.tot=cbind(sl.tot,mol.mRNA[,1],mol.scna[,1])

# Step II. Identifying survival-informative pairs using Cox regression
load("data/TCGA.drug.response/drug.patients.RData")		#drug.patients: TCGA samples that have treatment (drug) annocation
ix=which(!(prob$samples %in% drug.patients))			#we use the remainder of the samples to infer drug-SL network
numGenes = nrow(prob$scna[,ix])
numSamples = ncol(prob$scna[,ix])
mRNA.norm=prob$mRNA.norm[,ix]
scna.norm=prob$scna.norm[,ix]
surv.dt=prob$surv.dt[ix,]
age=qnorm.array(prob$age[ix])
sex=as.character(prob$sex[ix])
race=prob$race[ix]
types=prob$types[ix]
GII=qnorm.array(prob$GII[ix])
mRNAq2=prob$mRNAq2[,ix]
scnaq2=prob$scnaq2[,ix]

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
sl.tot=cbind(sl.tot,cox.sl.mRNA,cox.sl.scna)

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

sl.tot=cbind(sl.tot,pp,rank(pp,na.last="keep")/length(pp))
save(file="data/sl.pairs.patient.drug.response.RData",sl.tot)

