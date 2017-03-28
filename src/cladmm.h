//
// admm.h
// last update 01_12_2016
//

#ifndef __CLADMM_H
#define __CLADMM_H

#include <RcppArmadillo.h>

int cladmm(const arma::mat& input, arma::mat& z, arma::mat& u,
           const int& n, const double& admm_penalty, int maxiter, 
           const double& tolerance);

#endif