// Created by Zhifei Yan

#include "bfaddmin.h"
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

int bfaddmin(const mat& input, mat& z, mat& y, mat& u, mat& v, const int& n, 
             const int& k, const int& mmin, const double& admm_penalty, 
             int maxiter, const double& tolerance) {

  int niter;
	double rr, ss, tt, ww;
	uword num_neg;
	vec eigval(n);
	mat eigvec(n, n);
	mat x(n, n), y_old(n, n), z_old(n, n);
	vec temp(n), alpha(n);

	for(niter = 0; niter < maxiter; niter++) {
		// Store previous value of y and z
    y_old = y;
		z_old = z;

		// Projection onto the affine set
		x = 0.5 * (y + z - u - v + (input / admm_penalty));

		temp = 2.0 * sum(x, 1) - 2.0;
    
    alpha = 0.5 * temp / n;
    alpha -= 0.5 * ((trace(x) - k) / (n * (n - 1.0)) + 
                    sum(temp) * (n - 2.0) / (n * n * (2.0 * n - 2.0)));

    x.diag() -= 0.5 * (2.0 * (trace(x) - k) / (n - 1.0) - 
                       sum(temp) / (n * (n - 1.0)));
    x.each_col() -= alpha;
    x.each_row() -= alpha.t();
    
		// Projection onto the PSD cone
		y = x + v;
		eig_sym(eigval, eigvec, y);
		num_neg = sum(eigval < 0);

		// Reconstruct 
		if(num_neg == n) {
		  y.zeros();
		} else {
			y = (
			  eigvec.cols(num_neg, n - 1) * 
			  diagmat(eigval.subvec(num_neg, n - 1)) * 
			  eigvec.cols(num_neg, n - 1).t()
			);
    }

    // Projection onto the non-negative cone
    z = x + u;
    z.elem(find(z < 0)).zeros();
    z.elem(find(z > (1.0 / mmin))).fill(1.0 / mmin);
    
    // Dual variables update
    u = u + x - z;
    v = v + x - y;

    // Compute pieces constitute residual norms
    rr = norm(x - y, "fro");
    ss = norm(x - z, "fro");
    tt = norm(y - y_old, "fro");
    ww = norm(z - z_old, "fro");

    // Check convergence criterion
    if(rr < tolerance && ss < tolerance && tt < tolerance && ww < tolerance) {
      niter++;
      break;
    }
  }
  return niter;
}