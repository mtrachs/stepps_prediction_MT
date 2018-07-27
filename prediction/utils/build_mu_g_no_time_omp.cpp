// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_NO_DEBUG

#include <exception>
#include <iostream>
#include <RcppArmadillo.h>
#include <omp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export()]]
Rcpp::List build_mu_g_no_time_omp(NumericVector rho,
			  NumericVector mu_vec,
			  NumericVector alpha_s_vec,
			  NumericVector d_knots_vec,
			  NumericVector d_inter_vec,
			  int T, int K,
			  int od, int mu0,
			  NumericVector P_vec)
{
  IntegerVector dim_mu      = mu_vec.attr("dim");
  IntegerVector dim_alpha_s = alpha_s_vec.attr("dim");
  IntegerVector dim_d_knots = d_knots_vec.attr("dim");
  IntegerVector dim_d_inter = d_inter_vec.attr("dim");
  IntegerVector dim_P       = P_vec.attr("dim");

  mat mu(mu_vec.begin(), dim_mu[0], dim_mu[1], false);
  mat d_knots(d_knots_vec.begin(), dim_d_knots[0], dim_d_knots[1], false);
  mat d_inter(d_inter_vec.begin(), dim_d_inter[0], dim_d_inter[1], false);
  mat P(P_vec.begin(), dim_P[0], dim_P[1], false);

  int const W = K-1;
  int const niter = dim_mu[0];
  int const N = dim_d_inter[0];
  int const N_knots = dim_d_inter[1];

//   std::cout << "W " << W << std::endl;
  // std::cout << "N " << N << std::endl;
  // std::cout << "N_knots " << N_knots << std::endl;
   std::cout << "niter " << niter << std::endl;

  cube mu_g(N*T, W, niter);
  cube Halpha_s(N, W, niter);

  //#pragma omp parallel for
  for (int k=0; k<W; k++) {
    //for (int k=0; k<1; k++) {
    double rho_inv = 1.0 / rho[k];
    mat C_s = exp(-rho_inv * d_knots);
    mat c_s = exp(-rho_inv * d_inter);
    mat C_s_inv = inv(C_s);

    mat cs_Csinv;
    if (od) {
      cs_Csinv = c_s * C_s_inv - P * c_s * C_s_inv;
    } else {
      cs_Csinv = c_s * C_s_inv;
    }

    mat alpha_s(N_knots, niter);
    for (int i=0; i<niter; i++) {
      for (int v=0; v<N_knots; v++) {
        alpha_s(v, i) = alpha_s_vec[v*W*niter + k*niter + i];
      }
    }

    for (int i=0; i<niter; i++) {
      Halpha_s(0, k, i, size(N, 1, 1)) = cs_Csinv * alpha_s.col(i);

      if (mu0) {
	for (int j=0; j<N; j++) {
	  mu_g(j*T, k, i) = mu(i,k) + Halpha_s(j, k, i);
	}
      } else {
	throw std::runtime_error("mu0==0 not handled yet");
      }

      vec cs_Csinv_alpha_s = cs_Csinv * alpha_s.col(i);


      for (int t=1; t<T; t++) {

        if (mu0) {
	  for (int j=0; j<N; j++) {
	    mu_g(j*T+t, k, i) = mu(i,k) + cs_Csinv_alpha_s(j);
	  }
        } else {
	  throw std::runtime_error("mu0==0 not handled yet");
        }
      }

    }
  }

  return Rcpp::List::create(Rcpp::Named("mu_g") = mu_g,
			                      Rcpp::Named("mu") = mu,
			                      Rcpp::Named("Halpha_s") = Halpha_s);

}
