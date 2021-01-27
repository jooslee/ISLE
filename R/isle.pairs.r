# Identifying clinically relevant Synthetic LEthality (ISLE) for give gene pairs
# input: any candidate gene pair indices in N by 2 matrix format

isle.pairs <- function(sl){
# enable openMP support for gcc
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

load("data/prob.ach.new2.RData")
load("data/prob.ach.new1.RData")
load("data/prob.ach.old1.RData")
load("data/prob.mar.new1.RData")
load("data/prob.mar.old1.RData")
sourceCpp("src/essentialityTestPair.cpp",rebuild=T)
# col.names = (gene A SCNA, gene B ess) (gene A mRNA, gene B ess) 
# 			  (gene B SCNA, gene A ess) (gene B mRNA, gene A ess) for SL


# Step 0. define the initial pool
# (1) Initial Set I: Collect experimentally identified gold standard SL
load("data/sl.golden.set.RData")
gd$sr=gd$sl[gd$flag==1,]
sl0=gd$sr
ix=which(sl0[,1]>sl0[,2])
sl0[ix,]=sl0[ix,c(2,1)]

# (2) Initial Set II: Infer candidate SLs from single gene knockout data
rs=mclapply(seq(5),function(k){
	if (k==1) probc=prob.ach.new2
	if (k==2) probc=prob.ach.new1
	if (k==3) probc=prob.ach.old1		
	if (k==4) probc=prob.mar.new1
	if (k==5) probc=prob.mar.old1
	mRNAq=probc$mRNAq
	scnaq=probc$scnaq
	mRNAq[is.na(mRNAq)]=1
	scnaq[is.na(scnaq)]=1	
	mat=probc$mat
	p = essentialityTestPair(sl, t(scnaq), t(mRNAq), t(mat))
	return(p)
},mc.cores=32)
pcl=do.call(cbind,rs)
pcl[pcl< -999 | pcl==0.5] = NA
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
sl.tot=cbind(sl,p.ess,p.ess.adj)

# Step I. Identifying underrepresented pairs using hypergeometric test
load("data/prob.TCGA.RData")
sourceCpp("src/HyperGeometricTest.pair.cpp")
mol.scna = hypergeometricTestPair(scnaq=prob$scnaq2, pairs=sl)
mol.mRNA = hypergeometricTestPair(scnaq=prob$mRNAq2, pairs=sl)

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
sl.tot=cbind(sl.tot,mol.mRNA[,1],mol.scna[,1],cox.sl.mRNA,cox.sl.scna)

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
pp2=rank(pp,na.last="keep")/length(pp)

sl.tot=sl.tot=cbind(sl.tot,pp,pp2)
return(sl.tot)
}

