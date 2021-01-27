#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>		/* constants */
#define MIN(a,b) a<b?a:b 
#define MAX(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a
// #include <iostream>
// #include <fstream>



using namespace std;		// shorthand
using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec logRank(arma::mat times1, arma::mat times2){
	// cout << 25<< endl;
// events: 0 - death; 1 - censored
	int t1Len = times1.n_rows;
	int t2Len = times2.n_rows;
	// cout << 25<< endl;
	vec out(8);
	out.fill(-1000);
	if ((t1Len<1)|(t2Len <1))
	{
		return out;
	}
	// cout << 25<<d endl;
	mat times;
	double maxTime1 = max(times1.col(0));
	double maxTime2 = max(times2.col(0));
	times=  join_vert(times1, times2);
	uvec timesSort = sort_index(times.col(0));
// uvec timesUnique(t1Len + t2Len);
	// cout << 25<< endl;
	mat tab(t1Len + t2Len,9);
	tab.zeros();
	double timeq, timecurr, eventcurr;
	int timecnt;
 // 0:unique sorted time, 
 // 1: no of events in both groups, 
 // 2: no of events in  group 1, 
 // 3: no of events in  group 2, 
 // 4: no of death in  groups 1,  
 // 5: no of death in  groups 2,  
 // 6: no of people surviving in both groups,  
 // 7: no of people surviving in group 1, then fractional survival in group 1,  
 // 8: no of people surviving in group 2, then fractional survival in group 2,  
	timecnt =0; timeq = -10000;
	for (int tt = 0; tt < times.n_rows; tt++)
	{
		timecurr = times(timesSort(tt),0);
	// cout << 45<< endl;
		eventcurr = times(timesSort(tt),1);
		if (timeq != timecurr)
		{
			timeq = timecurr;
			timecnt++;
			tab(timecnt-1,0) = timeq; 
		}
		tab(timecnt-1,1) = tab(timecnt-1,1) +1; 
		if (eventcurr==0){
			tab(timecnt-1,5) = tab(timecnt-1,5) +1;  // total number of death that will change to death in group 2.... 
		} 
	// cout << 57<< endl;
	// cout << "tt" << tt << endl;
	} 
// cout << tab.rows(span(0,10)) << endl;
	// cout << 61<< endl;
	tab.resize(timecnt,9);

	// cout << 63<< endl;
	timesSort = sort_index(times1.col(0));

	timecnt =0; timeq = -10000;
	for (int tt = 0; tt < times1.n_rows; tt++)
	{
	// cout << 69<< endl;
		timecurr = times1(timesSort(tt),0);
		eventcurr = times1(timesSort(tt),1);
 	// cout << timecnt << " " << tab.n_rows << endl;
 	// cout << tab(0,0)<<endl;
		while (tab(timecnt,0) != timecurr)
		{
 	// cout << timecurr << endl;
 	// cout << timecnt << " " <<  tab(timecnt,0)<< endl;

			timecnt++; 
		}
 	// cout << timecnt << endl;
		tab(timecnt,2) = tab(timecnt,2) +1;
		if (eventcurr==0){
			tab(timecnt,4) = tab(timecnt,4) +1; 
		} 

	}

// cout << tab.rows(span(0,10)) << endl;
	tab.col(3) = tab.col(1) - tab.col(2);
	tab.col(5) = tab.col(5) - tab.col(4);
// cout << tab << endl;
	int tabSize = tab.n_rows;
	int npop;
// vec temp(tabSize); 

	npop = sum(tab.col(1));
	tab(0,6) = npop;
	if(tabSize < 2) return out;
	tab(span(1,tabSize-1), 6) = npop - cumsum(tab(span(0, tabSize -2),1));

	npop = sum(tab.col(2));
	tab(0,7) = npop;
	tab(span(1,tabSize-1), 7) = npop - cumsum(tab(span(0, tabSize -2),2));

	npop = sum(tab.col(3));
	tab(0,8) = npop;
	tab(span(1,tabSize-1), 8) = npop - cumsum(tab(span(0, tabSize -2),3));
// cout << tab << endl;
// tab = tab(tab(:,5) >0 | tab(:,6) >0,:);  % consider only deaths -> resize matrices to remove extras time when there were no deaths
// tab = tab(find((tab.col(4) >0) | (tab.col(5) >0),span(0,8)); 

	vec rate = (tab.col(4) + tab.col(5))/tab.col(6);
double Exp1 = sum( rate % tab.col(7));//% expected number of deaths group 1
double Exp2 = sum( rate % tab.col(8));//% expected number of deaths group 2


// double Exp1 = sum(tab.col(7));
// double Exp2 = sum(tab.col(8));
double Obs1 = sum(tab.col(4));
double Obs2 = sum(tab.col(5));


	// vec deltaTime1(tabSize), deltaTime2(tabSize);
double survivalRate1=1;
double survivalRate2=1;
double deltaTime1, deltaTime2;
double lastDeath1 = min(times1.col(0));
double lastDeath2 = min(times2.col(0));
	double startTime = 0; // min(lastDeath1, lastDeath2);
	lastDeath1 = lastDeath2  = startTime;

	double endTime = max(maxTime1, maxTime2);
	double auc1 =0, auc2 =0;
	for (int tt = 0; tt < tabSize; tt++){
		timecurr = tab(tt,0);
		if((tab(tt, 4) >0) || (timecurr == endTime) ){
			deltaTime1 = timecurr - lastDeath1;
			auc1 += deltaTime1 * survivalRate1;
			// deathRate1 = (tab(tt,4)/tab(tt,7)); 
			// deathRate1Abs = (tab(tt,4)/tab(tt,7))*survivalRate1; 
			survivalRate1 = survivalRate1 *( 1 -  (tab(tt,4)/tab(tt,7))); 
			lastDeath1 = timecurr; 
		}
		if((tab(tt, 5) >0) || (timecurr == endTime) ){
			deltaTime2 = timecurr - lastDeath2;
			auc2 += deltaTime2 * survivalRate2;
			survivalRate2 =  survivalRate2 *(1- (tab(tt,5)/tab(tt,8))); 
			lastDeath2 = timecurr; 
		}
		if(timecurr >=endTime)
			break;
	}
	auc1  = auc1/( endTime) ; 
	auc2  = auc2/( endTime) ; 

	//  = tab(span(1,tabSize -1), 0) -  tab(span(0,tabSize -2),0);
	// vec fractionalDeathT = tab.col(4)/tab.col(7);
	// vec fractionalDeath =  cumsum(fractionalDeathT); 
	// // double maxTime = tab(tabSize-1,0); // check tthis
	// double auc1 =  sum(fractionalDeath(span(0,tabSize-2)) % deltaTime)/maxTime;
	// auc1 = 1 -auc1;
	// fractionalDeathT = tab.col(5)/tab.col(8);
	// fractionalDeath =  cumsum(fractionalDeathT); 
	// double auc2 = sum(fractionalDeath(span(0,tabSize-2)) % deltaTime)/maxTime;
	// auc2 = 1 -auc2;

// We now calculate the significance of this result, using Chi distribution
	double chi_stat=((Exp1-Obs1)*(Exp1-Obs1))/Exp1+((Exp2-Obs2)*(Exp2-Obs2))/Exp2;
double p = 1- R::pchisq(chi_stat,1, 1, 0); //TODO

out(0) = p;
out(1) = Obs1;
out(2) = Exp1;
out(3) = Obs2;
out(4) = Exp2;
out(5) = chi_stat;
out(6)=auc1;
out(7)=auc2;

return out;
}





using namespace std;		// shorthand
using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat logRankMultiClass(arma::mat times, arma::uvec cls, int numClass){
	// cout << 25<< endl;
	// cls shout start from 1 
// events: 0 - death; 1 - censored
	int tLen = times.n_rows;
	// cout << 25<< endl;

	mat out(numClass, 5);
	out.fill(-1000);
	if (tLen<1)
	{
		return out;
	}
	// cout << 25<< endl;
	double maxTime = max(times.col(0));
	uvec timesSort = sort_index(times.col(0));
	// cout << 25<< endl;
	// mat tab(tLen,3 * (numClass + 1));
	// tab.zeros();
	mat eventMat(tLen, numClass + 1);
	mat deathMat(tLen, numClass + 1);
	vec timePoint(tLen);
	eventMat.zeros(); deathMat.zeros(); timePoint.zeros();
	double timeq, timecurr, eventcurr;
	int timecnt;
	int clscurr;
	// cout << trans(eventMat(span(0,10), 0)) << endl; 
	// cout << trans(deathMat(span(0,10), 0)) << endl; 
	
	timecnt =0; timeq = -100000;
	for (int tt = 0; tt < times.n_rows; tt++)
	{
		timecurr = times(timesSort(tt),0);
	// cout << 45<< endl;
		eventcurr = times(timesSort(tt),1);
		clscurr = cls(timesSort(tt)); /// this should start with because 0 corresponds to all groups
		if ((tt ==0 )|| (timeq != timecurr))
		{
			timeq = timecurr;
			timecnt++;
			timePoint(timecnt-1) = timeq; 
		}
		eventMat(timecnt-1,0) = eventMat(timecnt-1,0) +1; 
		eventMat(timecnt-1,clscurr) =  eventMat(timecnt-1,clscurr) +1; 
		if (eventcurr==0){
			deathMat(timecnt-1,0) =  deathMat(timecnt-1,0) +1; 
			deathMat(timecnt-1,clscurr) =  deathMat(timecnt-1,clscurr) +1; 
		}
	// cout << 57<< endl;
	// cout << "tt" << tt << endl;
	} 
// cout << tab.rows(span(0,10)) << endl;
	// cout << 61<< endl;
	// cout << "event " << trans(eventMat(span(0,10), 0)) << endl; 
	// cout << "death " << trans(deathMat(span(0,10), 0)) << endl; 
	// cout << "timePoint " << trans(timePoint(span(0,10))) << endl; 
	eventMat.resize(timecnt,numClass +1);
	deathMat.resize(timecnt,numClass +1);
	timePoint.resize(timecnt);
	mat survivalMat(timecnt, numClass + 1);
	vec chi_stat(numClass), p(numClass);
	// cout << 63<< endl;

	int tabSize = eventMat.n_rows;
	if(tabSize < 2) return out;
	int npop;
	npop = sum(eventMat.col(0));
	survivalMat(0,0) = npop;
	survivalMat(span(1,tabSize-1), 0) = npop - cumsum(eventMat(span(0, tabSize -2),0));
	vec rate = deathMat.col(0)/survivalMat.col(0);
	// cout <<  "survival" << trans(survivalMat(span(0,10), 0)) << endl; 
	// cout << "rate" << trans(rate(span(0,10))) << endl; 
	vec Exp(numClass), Obs(numClass);
	mat temp(1,1);
	for (int cc = 1; cc < numClass +1; cc++)
	{
		npop = sum(eventMat.col(cc));
		survivalMat(0,cc) = npop;
		survivalMat(span(1,tabSize-1), cc) = npop - cumsum(eventMat(span(0, tabSize -2),cc));
		 temp = trans(rate )* survivalMat.col(cc);//% expected number of deaths group 1
		 Exp(cc-1) =  temp(0,0);
		// cout << temp<< endl;
		 Obs(cc-1) = sum(deathMat.col(cc));
// We now calculate the significance of this result, using Chi distribution
		 chi_stat(cc-1)=((Exp(cc-1)-Obs(cc-1))*(Exp(cc-1)-Obs(cc-1)))/Exp(cc-1); 
		 p(cc-1) = 1- R::pchisq(chi_stat(cc-1),1, 1, 0); 
		}

		 // cout << temp.n_rows << temp.n_cols<< endl;

		vec survivalRate(numClass), deltaTime(numClass), auc(numClass);
		survivalRate.ones(); deltaTime.zeros();auc.zeros();
	// double startTime = 0; // min(lastDeath1, lastDeath2);
	// lastDeath1 = lastDeath2  = startTime;
		vec lastDeath(numClass); lastDeath.fill(timePoint(0) );

		double endTime = timePoint(tabSize -1);
		for (int tt = 0; tt < tabSize; tt++){
			timecurr = timePoint(tt);
		// cout << "tt" << tt << endl;
			for (int cc = 1; cc < numClass +1; cc++){

		// cout << cc << "cc" << tt << endl;
				if((deathMat(tt,cc) >0) || (timecurr == endTime) ){
					deltaTime(cc-1) = timecurr - lastDeath(cc-1);
					auc(cc-1) += deltaTime(cc-1) * survivalRate(cc-1);
					survivalRate(cc -1) = survivalRate(cc-1) *( 1 -  (deathMat(tt,cc)/survivalMat(tt,cc))); 
					lastDeath(cc-1) = timecurr; 
				}
			}
			if(timecurr >=endTime)
				break;
		}
		auc  = auc/( endTime) ; 
// out.col(0) = p;

		out.col(0) = p;
		out.col(1) = Obs;
		out.col(2) = Exp;
		out.col(3) = chi_stat;
		out.col(4)=auc;
		uvec nan_inx=  find(p != p );
		if(nan_inx.size() >0){
		// mat init(nan_inx.size(),1);
		// init.fill(-1000);
		// p(nan_inx).fill( -1000);
			out.rows(nan_inx).fill(-1000);
		}

		return out;
	}



using namespace std;		// shorthand
using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat logRankPairs(arma::umat pairs, arma::umat scna1Mat, arma::umat scna2Mat, arma::mat survival, arma::vec grps1, arma::vec grps2, int threads=1){
#ifdef _OPENMP
	if ( threads > 0 )
		omp_set_num_threads( threads ); 
#endif
	// cout << 145 << endl;
	int pairNum = pairs.n_rows;
	// cout << 145 << endl;
	int numSamples = scna1Mat.n_cols;
	// cout << 145 << endl;
	mat out1(pairNum, 8); 
	// cout << 145 << endl;
	int grps1Len = grps1.size();
	int grps2Len = grps2.size();
  #pragma omp parallel for schedule(static)
	for (int idx = 0; idx < pairNum; idx++)
	{
// int gene1, gene2;
// rowvec scna1, scna2;
	// cout << 145 << endl;
		uvec vDel(numSamples);
	// cout << 145 << endl;
		uvec vNorm(numSamples);

	// cout << 145 << endl;
            // gene1 = pairs(idx, 1);
            // gene2 = pairs(idx, 2);
	// cout << 166 << endl;
		urowvec scna1; 
		scna1 = scna1Mat.row(pairs(idx,0)-1);
	// cout << 169 << pairs(idx,1)<< endl;
		urowvec scna2;
		scna2 = scna2Mat.row(pairs(idx,1)-1);
		urowvec inx = 3 * scna2 + scna1 + 1;
		vDel.zeros(); vNorm.zeros();
		for (int i = 0; i < numSamples; ++i)
		{
			for (int jj = 0; jj < grps1Len; jj++)
				if(inx(i) == grps1(jj)) vDel(i) = 1;
			for (int jj = 0; jj < grps2Len; jj++)
				if(inx(i) == grps2(jj)) vNorm(i) = 1;
		}
            // vDel= find((inx ==2)|(inx ==3));
            // vNorm = find(inx==1);
	// cout << 172 << find(vDel)<< endl;
		uvec findVDel = find(vDel);
		mat times1(findVDel.size(),2);
		times1=survival.rows(findVDel); 
	// cout << 172 << times1.t() << endl;
		uvec findVNorm = find(vNorm);
		mat times2(findVNorm.size(),2);
		times2=survival.rows(findVNorm); 
	// cout << 172 << findVNorm<< endl;
	// cout << times2 << endl;
		out1.row(idx)= trans(logRank(times1,times2));     
	// cout << 192 << endl;
	}
	return out1;
}








using namespace std;		// shorthand
using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat aggregateLogRankPairs(arma::umat pairs, arma::umat scna1Mat, arma::umat scna2Mat, arma::mat survival, arma::vec grps1, arma::vec grps2, arma::vec typeInx, int typeNum=1, int threads=1){
#ifdef _OPENMP
	if ( threads > 0 )
		omp_set_num_threads( threads ); 
	#endif

	// cout << 145 << endl;
	int pairNum = pairs.n_rows;
	// cout << 145 << endl;
	int numSamples = scna1Mat.n_cols;
	// cout << 145 << endl;
	mat out(pairNum, 8); 
	// cout << 145 << endl;
	int grps1Len = grps1.size();
	int grps2Len = grps2.size();
  #pragma omp parallel for schedule(static)
	for (int idx = 0; idx < pairNum; idx++)
	{
// int gene1, gene2;
// rowvec scna1, scna2;

	// cout << 145 << endl;
            // gene1 = pairs(idx, 1);
            // gene2 = pairs(idx, 2);
	// cout << 166 << endl;
		urowvec scna1; 
		

		scna1 = scna1Mat.row(pairs(idx,0)-1);
		//cout << "here " << pairs(idx,0)<< endl;
		urowvec scna2;
		scna2 = scna2Mat.row(pairs(idx,1)-1);
		mat out1(typeNum, 8);
		int nzCnt =0; 
	for (int ii = 0; ii < typeNum; ii++)
	{
		// cout << ii << endl;
		uvec sel = arma::find(typeInx  ==  ii);
		// cout << 224 << endl;
		// cout << sel.t() << endl;

		urowvec scna1CurrType = scna1.cols(sel);
		urowvec scna2CurrType =  scna2.cols(sel);
		mat survivalCurrType = survival.rows(sel);
		urowvec inx = 3 * scna2CurrType + scna1CurrType + 1;


	// cout << 145 << endl;
	int numSamplesCurr = sel.size();
		uvec vDel(numSamplesCurr);
		uvec vNorm(numSamplesCurr);
		vDel.zeros(); vNorm.zeros();
	// cout << 145 << endl;
		for (int i = 0; i < numSamplesCurr; ++i)
		{
			for (int jj = 0; jj < grps1Len; jj++)
				if(inx(i) == grps1(jj)) vDel(i) = 1;
			for (int jj = 0; jj < grps2Len; jj++)
				if(inx(i) == grps2(jj)) vNorm(i) = 1;
		}
            // vDel= find((inx ==2)|(inx ==3));
            // vNorm = find(inx==1);
	// cout << 172 << find(vDel)<< endl;
		uvec findVDel = find(vDel);
		mat times1(findVDel.size(),2);
		times1=survivalCurrType.rows(findVDel); 
	// cout << 172 << times1.t() << endl;
		uvec findVNorm = find(vNorm);
		mat times2(findVNorm.size(),2);
		times2=survivalCurrType.rows(findVNorm); 
	// cout << 172 << findVNorm<< endl;
	// cout << times2 << endl;
	
		vec temp = logRank(times1,times2);
		// cout << 230 << endl;
		if(temp(0) >=0) {
			nzCnt ++;
		}else{
			temp.zeros();
		}
		// cout << 230 << endl;
		// cout << temp.t() << endl;
		// cout << out.n_rows<< "here " << ii << endl;
		out1.row(ii) = temp.t();
		//cout << temp.t()<< endl;

	// cout << 192 << endl;
	}
	// cout << 293 << endl;
	double Obs1 = sum( out1.col(1)); 
	double Exp1 = sum( out1.col(2)); 
	double Obs2 = sum( out1.col(3)); 
	double Exp2 = sum( out1.col(4)); 
	double auc1 = sum( out1.col(6)); 
	double auc2 = sum( out1.col(7)); 
	double chi_stat=((Exp1-Obs1)*(Exp1-Obs1))/Exp1+((Exp2-Obs2)*(Exp2-Obs2))/Exp2;
	double p = 1- R::pchisq(chi_stat,1, 1, 0); 
	// cout << 302 << endl;
	rowvec Tout(8); 
	Tout(0) = p;
	Tout(1) = Obs1;
	Tout(2) = Exp1;
	Tout(3) = Obs2;
	Tout(4) = Exp2;
	Tout(5) = chi_stat;
	Tout(6) = auc1/nzCnt;
	Tout(7) = auc2/nzCnt;
	// cout << 312 << endl;
	out.row(idx) = Tout;

	}
	return out;
}


using namespace std;		// shorthand
using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List aggregateLogRankPairsAll(arma::umat pairs, arma::umat scna1Mat, arma::umat scna2Mat, arma::mat survival, arma::vec grps1, arma::vec grps2, arma::vec typeInx, int typeNum=1, int threads=1){
#ifdef _OPENMP
	if ( threads > 0 )
		omp_set_num_threads( threads ); 
	#endif
	// cout << 145 << endl;
	int pairNum = pairs.n_rows;
	// cout << 145 << endl;
	int numSamples = scna1Mat.n_cols;
	// cout << 145 << endl;
	mat out(pairNum, 8); 
	mat pvalue(pairNum, typeNum); 
	mat auc(pairNum, typeNum); 
	// cout << 145 << endl;
	int grps1Len = grps1.size();
	int grps2Len = grps2.size();
  #pragma omp parallel for schedule(static)
	for (int idx = 0; idx < pairNum; idx++)
	{
// int gene1, gene2;
// rowvec scna1, scna2;

	// cout << 145 << endl;
            // gene1 = pairs(idx, 1);
            // gene2 = pairs(idx, 2);
	// cout << 166 << endl;
		urowvec scna1; 
		

		scna1 = scna1Mat.row(pairs(idx,0)-1);
	// cout << 169 << pairs(idx,1)<< endl;
		urowvec scna2;
		scna2 = scna2Mat.row(pairs(idx,1)-1);
		mat out1(typeNum, 8);
		int nzCnt =0; 
	for (int ii = 0; ii < typeNum; ii++)
	{
		// cout << ii << endl;
		uvec sel = arma::find(typeInx  ==  ii);
		// cout << 224 << endl;
		// cout << sel.t() << endl;

		urowvec scna1CurrType = scna1.cols(sel);
		urowvec scna2CurrType =  scna2.cols(sel);
		mat survivalCurrType = survival.rows(sel);
		urowvec inx = 3 * scna2CurrType + scna1CurrType + 1;


	// cout << 145 << endl;
	int numSamplesCurr = sel.size();
		uvec vDel(numSamplesCurr);
		uvec vNorm(numSamplesCurr);
		vDel.zeros(); vNorm.zeros();
	// cout << 145 << endl;
		for (int i = 0; i < numSamplesCurr; ++i)
		{
			for (int jj = 0; jj < grps1Len; jj++)
				if(inx(i) == grps1(jj)) vDel(i) = 1;
			for (int jj = 0; jj < grps2Len; jj++)
				if(inx(i) == grps2(jj)) vNorm(i) = 1;
		}
            // vDel= find((inx ==2)|(inx ==3));
            // vNorm = find(inx==1);
	// cout << 172 << find(vDel)<< endl;
		uvec findVDel = find(vDel);
		mat times1(findVDel.size(),2);
		times1=survivalCurrType.rows(findVDel); 
	// cout << 172 << times1.t() << endl;
		uvec findVNorm = find(vNorm);
		mat times2(findVNorm.size(),2);
		times2=survivalCurrType.rows(findVNorm); 
	// cout << 172 << findVNorm<< endl;
	// cout << times2 << endl;
	
		vec temp = logRank(times1,times2);
		// cout << 230 << endl;
		pvalue(idx, ii) = temp(0);
		auc(idx,ii) = temp(6) - temp(7);
		if(temp(0) >=0) {
			nzCnt ++;
		}else{
			temp.zeros();
		}
		// cout << 230 << endl;
		// cout << temp.t() << endl;
		// cout << out.n_rows<< "here " << ii << endl;
		out1.row(ii) = temp.t();
		//cout << temp.t()<< endl;

	// cout << 192 << endl;
	}
	// cout << 293 << endl;
	double Obs1 = sum( out1.col(1)); 
	double Exp1 = sum( out1.col(2)); 
	double Obs2 = sum( out1.col(3)); 
	double Exp2 = sum( out1.col(4)); 
	double auc1 = sum( out1.col(6)); 
	double auc2 = sum( out1.col(7)); 
	double chi_stat=((Exp1-Obs1)*(Exp1-Obs1))/Exp1+((Exp2-Obs2)*(Exp2-Obs2))/Exp2;
	double p = 1- R::pchisq(chi_stat,1, 1, 0); 
	// cout << 302 << endl;
	rowvec Tout(8); 
	Tout(0) = p;
	Tout(1) = Obs1;
	Tout(2) = Exp1;
	Tout(3) = Obs2;
	Tout(4) = Exp2;
	Tout(5) = chi_stat;
	Tout(6) = auc1/nzCnt;
	Tout(7) = auc2/nzCnt;
	// cout << 312 << endl;
	out.row(idx) = Tout;

	}
	  Rcpp::List outAll;
	  outAll["pancancer"] = out;
	  outAll["pvalue"] = pvalue;
	  outAll["auc"] = auc;
	return outAll;
}



using namespace std;            // shorthand
using namespace arma;           // shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat aggregateLogRankPairsAllCell(arma::umat pairs, arma::umat scna1Mat, arma::umat scna2Mat, arma::mat survival, arma::vec typeInx, int typeNum=1, int threads=1, int numCell = 9){
#ifdef _OPENMP
        if ( threads > 0 )
                omp_set_num_threads( threads );
                #endif
	int pairNum = pairs.n_rows;
        int numSamples = scna1Mat.n_cols;
        mat out(pairNum, 18);
  #pragma omp parallel for schedule(static)
        for (int idx = 0; idx < pairNum; idx++)
        {

                urowvec scna1;
                scna1 = scna1Mat.row(pairs(idx,0)-1);
                urowvec scna2;
                scna2 = scna2Mat.row(pairs(idx,1)-1);
                        rowvec Tout(18);
                for (int jj = 0; jj < numCell; jj++){
                        mat out1(typeNum, 8);
                        int nzCnt =0;
                        for (int ii = 0; ii < typeNum; ii++)
                        {
                                uvec sel = arma::find(typeInx  ==  ii);

                                urowvec scna1CurrType = scna1.cols(sel);
                                urowvec scna2CurrType =  scna2.cols(sel);
                                mat survivalCurrType = survival.rows(sel);
                                urowvec inx = 3 * scna2CurrType + scna1CurrType + 1;
                                int numSamplesCurr = sel.size();
                                uvec vDel(numSamplesCurr);
                                vDel.zeros();
                                for (int i = 0; i < numSamplesCurr; ++i)
                                        if(inx(i) == (jj +1)) vDel(i) = 1;
                                uvec findVDel = find(vDel);
                                mat times1(findVDel.size(),2);
                                times1=survivalCurrType.rows(findVDel);
                                uvec findVNorm = find(1-vDel);
                                mat times2(findVNorm.size(),2);
                                times2=survivalCurrType.rows(findVNorm);

                                vec temp = logRank(times1,times2);
                                if(temp(0) >=0) {
                                        nzCnt ++;
                                }else{
                                        temp.zeros();
                                }
                                out1.row(ii) = temp.t();

                        }
                        double Obs1 = sum( out1.col(1));
                        double Exp1 = sum( out1.col(2));
                        double Obs2 = sum( out1.col(3));
                        double Exp2 = sum( out1.col(4));
                        double auc1 = sum( out1.col(6));
                        double auc2 = sum( out1.col(7));
                        double chi_stat=((Exp1-Obs1)*(Exp1-Obs1))/Exp1+((Exp2-Obs2)*(Exp2-Obs2))/Exp2;
                        double p = 1- R::pchisq(chi_stat,1, 1, 0);
                        Tout((jj ) *2  ) = p;
                        Tout((jj) *2 +1  ) = (auc1 - auc2)/nzCnt;
                }
                out.row(idx) = Tout;

        }
        return(out);
}

using namespace std;            // shorthand
using namespace arma;           // shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List aggregateLogRankPairsMultiClass(arma::umat pairs, arma::umat scna1Mat, arma::umat scna2Mat, arma::mat survival, arma::vec typeInx, int typeNum=1, int threads=1, int numClass = 9){
#ifdef _OPENMP
	if ( threads > 0 )
		omp_set_num_threads( threads );
                #endif
	int pairNum = pairs.n_rows;
	int numSamples = scna1Mat.n_cols;
	mat ExpMat(pairNum, numClass), ObsMat(pairNum, numClass), aucMat(pairNum, numClass), pMat(pairNum,numClass);
	vec p_agg(pairNum);
  #pragma omp parallel for schedule(static)
	for (int idx = 0; idx < pairNum; idx++)
	{

		urowvec scna1;
		scna1 = scna1Mat.row(pairs(idx,0)-1);
		urowvec scna2;
		scna2 = scna2Mat.row(pairs(idx,1)-1);
		vec Exp(numClass), Obs(numClass), auc(numClass), chi_stat(numClass), p(numClass);
		Exp.zeros(); Obs.zeros(); auc.zeros(); chi_stat.zeros();p.zeros();
		vec nzCnt(numClass);
		nzCnt.zeros();
		for (int ii = 0; ii < typeNum; ii++)
		{
			uvec sel = arma::find(typeInx  ==  ii);
			// cout <<"sel" << sel.t() << endl;
			if(sel.size()> 0 ){


			urowvec scna1CurrType = scna1.cols(sel);
			urowvec scna2CurrType =  scna2.cols(sel);
			mat survivalCurrType = survival.rows(sel);
			urowvec inx = 3 * scna2CurrType + scna1CurrType + 1;
			// cout << survivalCurrType.t() << endl;
			mat temp = logRankMultiClass(survivalCurrType, inx.t(), numClass);
			// cout<< idx << "ii"<< ii << endl;
			for (int cc = 0; cc < numClass; cc++){
				if(temp(cc,0) >=0) {
					nzCnt(cc) = nzCnt(cc) +1;
				}else{
					temp.row(cc).fill(0);
				}
			}
			Exp = Exp + temp.col(2);
			Obs = Obs + temp.col(1);
			auc = auc + temp.col(4);
			}

		}
		// cout << 772 << endl;
		double chi_stat_agg = 0;
		for (int cc = 0; cc < numClass; cc++){
			if(nzCnt(cc) > 0 ) {
                auc(cc)  = auc(cc)/nzCnt(cc); //AUC normalization
            chi_stat(cc)=((Exp(cc)-Obs(cc))*(Exp(cc)-Obs(cc)))/Exp(cc); 
            chi_stat_agg += chi_stat(cc);
            p(cc) = 1- R::pchisq(chi_stat(cc), 1, 1, 0);  /// degree of freedom 0
            if(Exp(cc) < Obs(cc)) p(cc) = -p(cc);
                }else{
                	p(cc) = -1000;
                	auc(cc) = -1000;
                }	
        }
            p_agg(idx) = 1- R::pchisq(chi_stat_agg, numClass -1, 1, 0);  /// degree of freedom 0
        // cout << 780 << endl;

        ExpMat.row(idx) = Exp.t();
        ObsMat.row(idx) = Obs.t();
        pMat.row(idx) = p.t();
        aucMat.row(idx)=auc.t();

    }
    Rcpp::List out; 
    out["Exp"] = ExpMat; out["Obs"] = ObsMat; out["p"] = pMat; out["auc"] = aucMat;
    out["p_agg"] = p_agg;
    return(out);
}

using namespace std;            // shorthand
using namespace arma;           // shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List aggregateLogRankPairsSamples(arma::mat times1, arma::mat times2, arma::vec typeInx1, arma::vec typeInx2, int typeNum=1, int threads=1){
#ifdef _OPENMP
	if ( threads > 0 ) 
		omp_set_num_threads( threads ); 
               #endif
	vec out(8); 
	vec pvalue(typeNum); 
	vec auc(typeNum); 
	mat out1(typeNum, 8);
	int nzCnt =0;
  #pragma omp parallel for schedule(static)	
	for (int ii = 0; ii < typeNum; ii++)
	{
		uvec sel1 = arma::find(typeInx1  ==  ii);
		mat times1CurrType = times1.rows(sel1);
		uvec sel2 = arma::find(typeInx2  ==  ii);
		mat times2CurrType = times2.rows(sel2);		
		if ((times1CurrType.n_rows>1) && (times1CurrType.n_rows>1)){
		vec temp = logRank(times1CurrType,times2CurrType);
		pvalue(ii) = temp(0);
		auc(ii) = temp(6) - temp(7);
		if(temp(0) >=0) {
			nzCnt ++;
		}else{
			temp.zeros();
		}
		out1.row(ii) = temp.t();
	}
	}
	double Obs1 = sum( out1.col(1)); 
	double Exp1 = sum( out1.col(2)); 
	double Obs2 = sum( out1.col(3)); 
	double Exp2 = sum( out1.col(4)); 
	double auc1 = sum( out1.col(6)); 
	double auc2 = sum( out1.col(7)); 
	double chi_stat=((Exp1-Obs1)*(Exp1-Obs1))/Exp1+((Exp2-Obs2)*(Exp2-Obs2))/Exp2;
	double p = 1- R::pchisq(chi_stat,1, 1, 0); 
	vec Tout(8); 
	Tout(0) = p;
	Tout(1) = Obs1;
	Tout(2) = Exp1;
	Tout(3) = Obs2;
	Tout(4) = Exp2;
	Tout(5) = chi_stat;
	Tout(6) = auc1/nzCnt;
	Tout(7) = auc2/nzCnt;
	  Rcpp::List outAll;
	  outAll["pancancer"] = Tout;
	  outAll["pvalue"] = pvalue;
	  outAll["auc"] = auc;
	  outAll["out"] = out1;
	return outAll;
}
