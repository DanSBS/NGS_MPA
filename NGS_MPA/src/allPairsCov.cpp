#include <Rcpp.h>
#include <C:\cpplibs\armadillo-4.320.0\include\armadillo>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <omp.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP allPairsCov(SEXP R_x){
	omp_set_num_threads(8);
	NumericMatrix Rcpp_x(R_x);
	mat x(Rcpp_x.begin(), Rcpp_x.nrow(), Rcpp_x.ncol(), false);

	int n, nc;
	n = x.n_rows;
	nc = x.n_cols;
	mat o = ones<mat>(n, 1);
	mat cxy = zeros<mat>(nc, nc);

#pragma omp parallel
	{
#pragma omp for
		for(int i=0; i<n; i++){
			mat xrep = o * x.row(i);
			cxy += xrep.t() * x;
		}
	}

	return(Rcpp::wrap(cxy));

}
