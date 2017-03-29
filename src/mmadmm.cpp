//
// mmadmm.cpp
// created by Zhifei Yan
//

#include "mmadmm.h"
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

int mmadmm(const mat& input, mat& z, mat& u, const int& n, 
           const double& admm_penalty, int maxiter, 
           const double& tolerance) {
  int niter;
  double rr, ss;
  uword num_neg;
  vec eigval(n);
  mat eigvec(n, n);
  mat y(n, n), z_old(n, n);

  for(niter = 0; niter < maxiter; niter++) {
    // Store previous value of z
    z_old = z;

    // Projection onto the PSD cone
    y = z - u - input / admm_penalty;
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

    // Projection onto the elementwise 0 to 1 bound and diagonal entries be 1
    z = y + u;
    z.elem(find(z < 0)).zeros();
    z.elem(find(z > 1)).ones();
    z.diag().ones();

    // Dual variables update
    u = u + y - z;

    // Compute pieces constitute residual norms
    rr = norm(z - y, "fro");
    ss = norm(z - z_old, "fro");

    // Check convergence criterion
    if(rr < tolerance && ss < tolerance) {
      niter++;
      break;
    }
  }

  return niter;
}