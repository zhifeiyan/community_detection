#include "gvadmm.h"
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

int gvadmm(const mat& input, mat& z, mat& y, mat& u, mat& v, const int& n, 
         const double& lambda, const double& admm_penalty, int maxiter, 
         const double& tolerance) {

  int niter;
	double rr, ss, tt, ww;
	uword num_neg;
	vec eigval(n);
	mat eigvec(n, n);

	mat x(n, n), y_old(n, n), z_old(n, n);
	
	double sum_diag, sum_all;

	for(niter = 0; niter < maxiter; niter++) {
		// Store previous value of y and z
    y_old = y;
		z_old = z;

		// Projection onto the affine set
		x = 0.5 * (y + z - u - v + (input / admm_penalty));

    sum_diag = sum(x.diag());
    sum_all = accu(x);

    x.diag().fill((sum_all + std::pow(n, 2) - sum_diag - lambda) / 
                  (std::pow(n, 2) - n));
    x -= (sum_diag - n - sum_all + lambda) / (n - std::pow(n, 2));

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