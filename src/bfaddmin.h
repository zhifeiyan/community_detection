// Created by Zhifei Yan

#ifndef __BFADDMIN_H
#define __BFADDMIN_H

#include <RcppArmadillo.h>

int bfaddmin(const arma::mat& input, arma::mat& z, arma::mat& y, 
             arma::mat& u, arma::mat& v, const int& n, const int& k,
             const int& mmin,  const double& admm_penalty, 
             int maxiter, const double& tolerance);

#endif