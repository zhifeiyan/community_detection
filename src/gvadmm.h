//
// gvadmm.cpp
// created by Zhifei Yan
//

#ifndef __GVADMM_H
#define __GVADMM_H

#include <RcppArmadillo.h>

int gvadmm(const arma::mat& input, arma::mat& z, arma::mat& y, 
           arma::mat& u, arma::mat& v, const int& n, const double& lambda, 
           const double& admm_penalty, int maxiter, const double& tolerance);

#endif