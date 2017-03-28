#include <RcppArmadillo.h>
#include <cmath>

#include "gvadmm.h"

using namespace Rcpp;
using namespace arma;

//' SDP-GV method
//'
//' This function computes a solution path of SDP-GV
//' under a sequence of values of tuning paramter lambda
//' (see COMMUNITY DETECTION IN SPARSE NETWORKS VIA GROTHENDIECK’S INEQUALITY by OLIVIER GUE ́DON AND ROMAN VERSHYNIN)
//'
//' @param S Input symmetric matrix
//' @param lambda Vector of possible values of SDP tuning parameter
//' @param maxiter Max number of iteration
//' @param tolerance ADMM tolerance
//' @param admm_penalty ADMM penalty parameter
//' @param verbose Level of verbosity
//' @return A list containing SDP results
//' @export
// [[Rcpp::export]]
List sdp_gv(NumericMatrix S, NumericVector lambda = NumericVector::create(),
            int maxiter = 1e6, double tolerance = 1e-2,
            double admm_penalty = 100.0, int verbose = 0) {

	// Sanity checks
	if(S.nrow() < 2) stop("Expected S to be a matrix");
	if(maxiter < 1) stop("Expected maxiter > 0");
	if(tolerance <= 0.0) stop("Expected tolerance > 0");

	int nsol;
	if(lambda.size() > 0) {
		nsol = lambda.size();
	} else {
		stop("Expected length of lambda > 0");
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

	// Outer loop to compute the solution path
	for(int i = 0; i < nsol; i++) {
		if(verbose > 0) Rcout << ".";

		// ADMM
		niter[i] = gvadmm(_S, z, y, u, v, ndim, lambda[i], 
                      admm_penalty, maxiter, tolerance);

		// Store solution
		cluster[i] = z;
		maxval[i] = dot(_S, z);

		if(verbose > 1) Rcout << niter[i];
	}

	if(verbose > 0) Rcout << std::endl;

	// Return
	List out = List::create(
	  Named("lambda") = lambda,
		Named("cluster") = cluster,
		Named("maxval") = maxval,
		Named("niter") = niter
	);

  return out;
}
