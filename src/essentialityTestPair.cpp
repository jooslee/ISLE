#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>      /* constants */
#define MIN(a,b) a<b?a:b 
#define MAX(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a
// #include <iostream>
// #include <fstream>

using namespace std;        // shorthand
using namespace arma;       // shorthand
using namespace Rcpp;       // shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double  myRanksum(arma::vec x, arma::vec y) {
    //alternative hypothesis x > y 
 // x and y are assumed to be sorted. 
 int xlen = x.size();
 int ylen = y.size();
 if (!((xlen >0) & (ylen > 0)))
    return -1000;
 vec xy = join_vert(x,y);
 uvec sortxy = sort_index(xy); 
 // arma uvec yi = sort_index(y); 
 // arma vec rankx(xlen);
 // arma vec ranky(ylen);
 double ranksumx = 0, ranksumy =0;
 double curr, next;
 int currCount =-1;
 double xinacc =0;
 vec temp;
 double xcnt=0; double ycnt =0;
 double acc =0;
 double pval;
 for (int i = 0; i < xlen + ylen; ++i)
 {
    curr = xy(sortxy(i));
    if(i < ((xlen + ylen) -1) )
        next = xy(sortxy(i +1));
    else 
        next = -1000;

    currCount++;
    if(sortxy(i) < xlen) 
        xcnt++;
    else 
        ycnt++;

    acc = acc + i + 1;
    
    if(curr != next){
        xinacc =   acc * xcnt/(xcnt + ycnt);   
        ranksumx = ranksumx + xinacc;
        // if(ranksumx != ranksumx){
            // cout<< ranksumx << " " << xinacc << " " << acc << endl;
            // return 0;
        // }
        // cout << ranksumx<<endl;
        acc = 0;
        xcnt =0; ycnt =0;
    }
    // cout <<xinacc << " " <<  sortxy(i) << " " << acc << " " << ranksumx << endl;
     
 }
 // cout <<  sortxy.t() << endl;
// cout << ranksumx << endl;
// cout << xlen << endl;
// cout << xlen*(xlen +1)/2.0 << endl;
// if sample size is less than 10 use pwilcox otherwise use normal approximation
 // cout<< ranksumx << " " << xinacc << " " << acc << " "<< xlen << " " << ylen << endl;
 if ((xlen <= 9) & (ylen <=9))
 {
 ranksumx= ranksumx - xlen*(xlen +1)/2.0;
    pval = R::pwilcox(ranksumx,  xlen,  ylen, 0, 0);
     /* code */
 } else{
    double mu =  xlen* (xlen + ylen +1) /2.0; 
    // double sigma = sqrt(xlen * ylen * (xlen  + ylen +1 )/12.0 ); /// this causes overflow of numbers
    double sigma =  sqrt(ylen) * sqrt(mu/6.0);
    // cout << mu << " " << sigma << endl; 
    pval = R::pnorm(ranksumx,  mu,  sigma, 0, 0);

 }
 return pval;
}

using namespace Rcpp;       // shorthand
using namespace arma;       // shorthand
using namespace std;        // shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat essentialityTest(arma::uvec scna1, arma::uvec mRNA1, arma::mat scnaq, arma::mat mRNAq, arma::mat ess, int threads=1){
   if ( threads > 0 )
    omp_set_num_threads( threads ); 
    int numGenes = scnaq.n_cols;
    mat out(numGenes,2);
     // scna2, mRNA2;
    // umat cnt(3,3);
 //cout << 102 << endl;
    uvec grp1, grp2, grp3, grp4;
    grp1 = find(scna1 ==0);
     //cout << grp1.t() << endl;
    grp2 = find(mRNA1 ==0);
    grp3 = find(scna1 ==2);
    grp4 = find(mRNA1 ==2);
    int gene2;
  // #pragma omp parallel for 
  #pragma omp parallel for schedule(static)
    for ( gene2 = 0; gene2< numGenes; gene2++)
    {
 //cout << 106 << endl;
     ////   vec scna2 = scnaq.col(gene2);
 //cout << 108 << " " << scna2.size()<< endl;
    //cout << "gene2" << gene2 << "\n" << endl;
        vec ess1 = ess.col(gene2);
        rowvec temp(2);
        temp( 0) = myRanksum(ess1(grp3), ess1(grp1));
	//cout << "1" << endl
        temp( 1) = myRanksum(ess1(grp4), ess1(grp2));
	//cout << "2" << endl
       //// temp( 2) = myRanksum(scna2(grp2), scna2(grp3));
	//cout << "3" << endl
////        temp( 3) = myRanksum(scna2(grp1), scna2(grp4));
	//cout << "4" << endl
    ////    temp( 4) = myRanksum(mRNA2(grp1), mRNA2(grp2));
	//cout << "5" << endl
 ////       temp( 5) = myRanksum(mRNA2(grp1), mRNA2(grp3));
	//cout << "6" << endl
     ////   temp( 6) = myRanksum(mRNA2(grp2), mRNA2(grp3));
	//cout << "7" << endl
////        temp( 7) = myRanksum(mRNA2(grp1), mRNA2(grp4));
	//cout << "8" << endl
        // temp( 7) = gene2;
        // myRanksum(mRNA2(grp1), mRNA2(grp4));
        // #pragma omp critical
        // {
        //cout << temp << endl;
            out.row(gene2 ) = temp; 
        // }
    }
    return out;
}

using namespace Rcpp;       // shorthand
using namespace arma;       // shorthand
using namespace std;        // shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat essentialityTestAll(arma::mat scnaq, arma::mat mRNAq, arma::mat ess, int threads=1){
   if ( threads > 0 )
    omp_set_num_threads( threads ); 
    int numGenes = scnaq.n_cols;
    mat out(numGenes*numGenes,2);
    int gene1, gene2;
    uvec grp1, grp2, grp3, grp4;
    for ( gene1 = 0; gene1 < numGenes; gene1++ )
    {
    vec mRNA1=mRNAq.col(gene1);
    vec scna1=scnaq.col(gene1);
    grp1 = find(scna1 ==0);
    grp2 = find(mRNA1 ==0);
    grp3 = find(scna1 ==2);
    grp4 = find(mRNA1 ==2);
  #pragma omp parallel for schedule(static)
    for ( gene2 = 0; gene2< numGenes; gene2++) if (gene1!=gene2)
    {
        vec ess1 = ess.col(gene2);
        rowvec temp(2);
        temp( 0) = myRanksum(ess1(grp3), ess1(grp1));
        temp( 1) = myRanksum(ess1(grp4), ess1(grp2));
            out.row(numGenes*(gene1)+gene2 ) = temp; 
    }}
    return out;
}



using namespace arma;       // shorthand
using namespace Rcpp;		// shorthand
using namespace std;		// shorthand
// [[Rcpp::export]]
arma::mat essentialityTestPair(Rcpp::IntegerMatrix pairs, arma::mat scnaq, arma::mat mRNAq, arma::mat ess, int threads=1){
	int numPairs = pairs.nrow();
	int gene1, gene2;
	
	int xx;
    mat out(numPairs,4);

 if ( threads > 0 )
    omp_set_num_threads( threads );
	for (int pairInx = 0; pairInx< numPairs; pairInx++)
	{
		gene1 = pairs(pairInx,0)-1;
		gene2 = pairs(pairInx,1)-1;

	    uvec grp1, grp2, grp3, grp4, grq1, grq2, grq3, grq4;
	
	    vec mRNA1=mRNAq.col(gene1);	    
	    vec scna1=scnaq.col(gene1);
	    vec mRNA2=mRNAq.col(gene2);	    
	    vec scna2=scnaq.col(gene2);

	    grp1 = find(scna1 ==0);	    grp2 = find(mRNA1 ==0);
	    grp3 = find(scna1 ==2);	    grp4 = find(mRNA1 ==2);
	    grq1 = find(scna2 ==0);	    grq2 = find(mRNA2 ==0);
	    grq3 = find(scna2 ==2);	    grq4 = find(mRNA2 ==2);
	    
        vec ess1 = ess.col(gene2);
        vec ess2 = ess.col(gene1);        
        rowvec temp(4);
        temp( 0) = myRanksum(ess1(grp3), ess1(grp1));
        temp( 1) = myRanksum(ess1(grp4), ess1(grp2));
        temp( 2) = myRanksum(ess2(grq3), ess2(grq1));
        temp( 3) = myRanksum(ess2(grq4), ess2(grq2));
            out.row(pairInx ) = temp; 
    }
    return out;
}
