//
// mmadmm.h
// created by Zhifei Yan
//

#ifndef __MMADMM_H
#define __MMADMM_H

#include <RcppArmadillo.h>

int mmadmm(const arma::mat& input, arma::mat& z, arma::mat& u,
           const int& n, const double& admm_penalty, int maxiter, 
           const double& tolerance);

#endif