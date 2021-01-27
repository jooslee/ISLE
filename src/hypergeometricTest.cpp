#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>		/* constants */
#define MIN(a,b) a<b?a:b 
#define MAX(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a


using namespace Rcpp;		// shorthand
using namespace std;		// shorthand

Rcpp::NumericMatrix hypergeometricTest(Rcpp::IntegerVector scna1, Rcpp::IntegerMatrix scnaq, int lowerTail=1){
	int numGenes = scnaq.nrow();
	int numSamples = scnaq.ncol();

	int numSamplesnz =0;
	NumericMatrix hyperpscna(numGenes,9);
	IntegerVector scna2;
	IntegerVector counts(9);

	IntegerVector inx(numSamples);
	IntegerVector gene1cnt(3), gene2cnt(3);
	IntegerVector genecnt(9);
	NumericVector p(9);
	int xx;

	for (int gene2 = 0; gene2< numGenes; gene2++)
	{
		scna2 = scnaq(gene2,_);
		inx = 3 * scna2 + scna1;
		std::fill(counts.begin(), counts.end(), 0);
		std::fill(gene1cnt.begin(), gene1cnt.end(), 0);
		std::fill(gene2cnt.begin(), gene2cnt.end(), 0);
		inx = 3 * scna2 + scna1;

		numSamplesnz =0;	
		for (int samp = 0; samp < numSamples; ++samp)
		{

			if((scna1[samp] >=0) & (scna2[samp] >=0)) {
				numSamplesnz ++;
				counts[inx[samp]]++;
			}
		}	

		for (int row = 0; row < 3; row++)
		{
			for (int col = 0; col < 3; col++)
			{
				xx = row *3 + col;
				gene1cnt(col) = gene1cnt(col) + counts(xx);
				gene2cnt(row) = gene2cnt(row) + counts(xx);
			}
		}

	for (int row = 0; row < 3; row++)
		{
			for (int col = 0; col < 3; col++)
			{
				xx = row *3 + col;
			hyperpscna(gene2,xx) = R::phyper(counts(xx), gene1cnt(col), numSamplesnz - gene1cnt(col),  gene2cnt(row), lowerTail, 0);
			}
		}
	}
	return hyperpscna;

}