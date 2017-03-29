//
// sdp_caili.cpp
// created by Zhifei Yan
//

#include <RcppArmadillo.h>
#include <cmath>

#include "cladmm.h"

using namespace Rcpp;
using namespace arma;

//' SDP-CL method
//'
//' This function computes a solution path of SDP-CL 
//' under a sequence of values of tuning parameter lambda
//' (see ROBUST AND COMPUTATIONALLY FEASIBLE COMMUNITY DETECTION IN THE PRESENCE OF ARBITRARY OUTLIER NODES by T. Tony Cai and Xiaodong Li)
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
List sdp_caili(NumericMatrix S, NumericVector lambda = NumericVector::create(),
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

  // Wrap the input adjacency matrix with an arma::mat
  const mat _S(S.begin(), ndim, ndim, false);

  // Placeholders for solutions
  List cluster(nsol);
  IntegerVector niter(nsol);
  // NumericVector maxval(nsol);

  // ADMM variables, passing by reference to the ADMM algorithm
  // The latter case uses the results of the previous case as initial values
  // z keeps track of the solution matrix
  mat z = zeros<mat>(ndim, ndim),
      u = zeros<mat>(ndim, ndim),
      input;

  // Outer loop to compute the solution path
  for(int i = 0; i < nsol; i++) {
    if(verbose > 0) Rcout << ".";

    // Generate input matrix of SDP using current value of lambda
    // input matrix = - adj + lambda (E - I)
    input = -_S + lambda(i);
    input.diag() -= lambda(i);

    // ADMM
    niter[i] = cladmm(input, z, u, ndim, admm_penalty, maxiter, tolerance);

    // Store solution
    cluster[i] = z;
    // maxval[i] = dot(_S, z);

    if(verbose > 1) Rcout << niter[i];
  }

  if(verbose > 0) Rcout << std::endl;

  // Return
  List out = List::create(
    Named("lambda") = lambda,
    Named("cluster") = cluster,
    // Named("maxval") = maxval,
    Named("niter") = niter
  );

  return out;
}
