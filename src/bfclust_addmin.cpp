// Created by Zhifei Yan

#include <RcppArmadillo.h>
#include <cmath>

#include "bfaddmin.h"

using namespace Rcpp;
using namespace arma;

//' SDP-bf adding an entry upper bound constraint
//'
//' @param S Input symmetric matrix
//' @param mmin Minimum cluster size
//' @param nclust Vector of candidate numbers of clusters
//' @param maxiter Max number of iteration
//' @param tolerance ADMM tolerance
//' @param admm_penalty ADMM penalty parameter
//' @param verbose Level of verbosity
//' @return A list containing SDP results
//' @export
// [[Rcpp::export]]
List bfclust_addmin(NumericMatrix S, int mmin, 
                    IntegerVector nclust = IntegerVector::create(),
                    int maxiter = 1e2, double tolerance = 1e-2,
                    double admm_penalty = 100.0, int verbose = 0) {

	// Sanity checks
	if(S.nrow() < 2) stop("Expected S to be a matrix");
	if(maxiter < 1) stop("Expected maxiter > 0");
	if(tolerance <= 0.0) stop("Expected tolerance > 0");

	int nsol;
	if(nclust.size() > 0) {
		nsol = nclust.size();
	} else {
		stop("Expected length of nclust > 0");
	}

	int ndim = S.nrow();

	// Wrap the input matrix with an arma::mat
	const mat _S(S.begin(), ndim, ndim, false);

	// Placeholders for solutions
	List cluster(nsol);
	IntegerVector niter(nsol);
	NumericVector maxval(nsol);

	// ADMM variables, passing by reference to the ADMM algorithm
	// The latter case uses the results of the previous case as initial values
	// z keeps track of the solution matrix
	mat z = zeros<mat>(ndim, ndim),
      y = zeros<mat>(ndim, ndim),
		  u = zeros<mat>(ndim, ndim),
		  v = zeros<mat>(ndim, ndim);

	// ADMM parameters
	// double tolerance_abs = std::sqrt(ndim) * tolerance;

	// Outer loop to compute the solution path
	for(int i = 0; i < nsol; i++) {
		if(verbose > 0) Rcout << ".";

		// ADMM
		niter[i] = bfaddmin(_S, z, y, u, v, ndim, nclust[i], mmin, 
                        admm_penalty, maxiter, tolerance);

		// Store solution
		cluster[i] = z;
		maxval[i] = dot(_S, z);

		if(verbose > 1) Rcout << niter[i];
	}

	if(verbose > 0) Rcout << std::endl;

	// Return
	List out = List::create(
	  Named("nclust") = nclust,
		Named("cluster") = cluster,
		Named("maxval") = maxval,
		Named("niter") = niter
	);

  return out;
}
