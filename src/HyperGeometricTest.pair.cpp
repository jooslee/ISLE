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
// #include <iostream>
// #include <fstream>


using namespace Rcpp;		// shorthand
using namespace std;		// shorthand
// [[Rcpp::export]]
Rcpp::NumericMatrix hypergeometricTest(Rcpp::IntegerVector scna1, Rcpp::IntegerMatrix scnaq, int lowerTail=1){
	int numGenes = scnaq.nrow();
	int numSamples = scnaq.ncol();
	// double numSamples2 =numSamples * numSamples;
	int numSamplesnz =0;
	NumericMatrix hyperpscna(numGenes,9);
	IntegerVector scna2;
	IntegerVector counts(9);
	// IntegerMatrix cnt(3,3);
	IntegerVector inx(numSamples);
	IntegerVector gene1cnt(3), gene2cnt(3);
	IntegerVector genecnt(9);
	NumericVector p(9);
	int xx;
	// cout << numSamples<< endl;
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
			// cout <<  inx[samp] << endl;
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


using namespace Rcpp;		// shorthand
using namespace std;		// shorthand
// [[Rcpp::export]]
Rcpp::NumericMatrix hypergeometricTestPair(Rcpp::IntegerMatrix pairs, Rcpp::IntegerMatrix scnaq, int lowerTail=1, int threads=1){
	// if ( threads > 0 )
		// omp_set_num_threads( threads );	 
	int numPairs = pairs.nrow();
	int numSamples = scnaq.ncol();
	int numSamplesnz =0;
	NumericMatrix hyperpscna(numPairs,9);
	// double numSamples2;// =numSamples * numSamples;
	IntegerVector scna2;
	IntegerVector scna1;
	IntegerVector counts(9);
	IntegerVector inx(numSamples);
	IntegerVector gene1cnt(3), gene2cnt(3);
	IntegerVector genecnt(9);
	IntegerVector scnaInx;
	NumericVector p(9);
	cout << "here" << endl;
	int xx;
  // #pragma omp parallel for schedule(static)
	for (int pairInx = 0; pairInx< numPairs; pairInx++)
	{
		scna1 = scnaq(pairs(pairInx,0)-1,_);
		scna2 = scnaq(pairs(pairInx,1)-1,_);
		std::fill(counts.begin(), counts.end(), 0);
		std::fill(gene1cnt.begin(), gene1cnt.end(), 0);
		std::fill(gene2cnt.begin(), gene2cnt.end(), 0);
		// scnaInx = find((scna1 >=0) | (scna2 >=0)  );
		// scna1 = scna1(scnaInx);
		// scna1 = scna1(scnaInx);
		inx = 3 * scna2 + scna1;
		std::fill(counts.begin(), counts.end(), 0);
		std::fill(gene1cnt.begin(), gene1cnt.end(), 0);
		std::fill(gene2cnt.begin(), gene2cnt.end(), 0);
		numSamplesnz =0;
		// cout << inx << endl;

		for (int samp = 0; samp < numSamples; ++samp)
		{
			// cout <<  inx[samp] << endl;
			if((scna1[samp] >=0) & (scna2[samp] >=0)) {
				numSamplesnz ++;
				counts[inx[samp]]++;
			}
		}

		// numSamples2 = numSamplesnz * numSamplesnz;

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
			hyperpscna(pairInx,xx) = R::phyper(counts(xx), gene1cnt(col), numSamplesnz - gene1cnt(col),  gene2cnt(row), lowerTail, 0);
			}
		}
	}
	return hyperpscna;

}


using namespace Rcpp;		// shorthand
using namespace std;		// shorthand
// [[Rcpp::export]]
Rcpp::NumericMatrix hypergeometricTestPairmiRNA(Rcpp::IntegerMatrix pairs, Rcpp::IntegerMatrix scnaq, Rcpp::IntegerMatrix miRNAq, int lowerTail=1, int threads=1){
	// if ( threads > 0 )
		// omp_set_num_threads( threads );	 
	int numPairs = pairs.nrow();
	int numSamples = scnaq.ncol();
	int numSamplesnz =0;
	NumericMatrix hyperpscna(numPairs,9);
	// double numSamples2;// =numSamples * numSamples;
	IntegerVector scna2;
	IntegerVector scna1;
	IntegerVector counts(9);
	IntegerVector inx(numSamples);
	IntegerVector gene1cnt(3), gene2cnt(3);
	IntegerVector genecnt(9);
	IntegerVector scnaInx;
	NumericVector p(9);
	cout << "here" << endl;
	int xx;
  // #pragma omp parallel for schedule(static)
	for (int pairInx = 0; pairInx< numPairs; pairInx++)
	{
		scna1 = scnaq(pairs(pairInx,0)-1,_);
		scna2 = miRNAq(pairs(pairInx,1)-1,_);
		std::fill(counts.begin(), counts.end(), 0);
		std::fill(gene1cnt.begin(), gene1cnt.end(), 0);
		std::fill(gene2cnt.begin(), gene2cnt.end(), 0);
		// scnaInx = find((scna1 >=0) | (scna2 >=0)  );
		// scna1 = scna1(scnaInx);
		// scna1 = scna1(scnaInx);
		inx = 3 * scna2 + scna1;
		std::fill(counts.begin(), counts.end(), 0);
		std::fill(gene1cnt.begin(), gene1cnt.end(), 0);
		std::fill(gene2cnt.begin(), gene2cnt.end(), 0);
		numSamplesnz =0;
		// cout << inx << endl;

		for (int samp = 0; samp < numSamples; ++samp)
		{
			// cout <<  inx[samp] << endl;
			if((scna1[samp] >=0) & (scna2[samp] >=0)) {
				numSamplesnz ++;
				counts[inx[samp]]++;
			}
		}

		// numSamples2 = numSamplesnz * numSamplesnz;

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
			hyperpscna(pairInx,xx) = R::phyper(counts(xx), gene1cnt(col), numSamplesnz - gene1cnt(col),  gene2cnt(row), lowerTail, 0);
			}
		}
	}
	return hyperpscna;

}

