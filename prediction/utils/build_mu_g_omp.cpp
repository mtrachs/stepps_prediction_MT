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
Rcpp::List build_mu_g_omp(NumericVector rho,
			  NumericVector sigma_vec,
			  NumericVector lambda_vec,
			  NumericVector mu_vec,
			  NumericVector mu_t_vec,
			  NumericVector alpha_s_vec,
			  NumericVector alpha_t_vec,
			  NumericVector d_knots_vec,
			  NumericVector d_inter_vec,
			  int T, int K,
			  int od, int mu0,
			  NumericVector P_vec)
{
  IntegerVector dim_sigma   = sigma_vec.attr("dim");
  IntegerVector dim_lambda  = lambda_vec.attr("dim");
  IntegerVector dim_mu      = mu_vec.attr("dim");
  IntegerVector dim_mu_t    = mu_t_vec.attr("dim");
  IntegerVector dim_alpha_s = alpha_s_vec.attr("dim");
  IntegerVector dim_alpha_t = alpha_t_vec.attr("dim");
  IntegerVector dim_d_knots = d_knots_vec.attr("dim");
  IntegerVector dim_d_inter = d_inter_vec.attr("dim");
  IntegerVector dim_P       = P_vec.attr("dim");

  mat mu(mu_vec.begin(), dim_mu[0], dim_mu[1], false);
  mat lambda(lambda_vec.begin(), dim_lambda[0], dim_lambda[1], false);
  mat d_knots(d_knots_vec.begin(), dim_d_knots[0], dim_d_knots[1], false);
  mat d_inter(d_inter_vec.begin(), dim_d_inter[0], dim_d_inter[1], false);
  mat P(P_vec.begin(), dim_P[0], dim_P[1], false);

  int const W = K-1;
  int const niter = dim_sigma[0];
  int const N = dim_d_inter[0];
  int const N_knots = dim_d_inter[1];

  // std::cout << "W " << W << std::endl;
  // std::cout << "N " << N << std::endl;
  // std::cout << "N_knots " << N_knots << std::endl;
  // std::cout << "niter " << niter << std::endl;

  cube mu_g(N*T, W, niter);
  cube Halpha_s(N, W, niter);
  cube Halpha_t;

  if (mu0) {
    Halpha_t.resize(N, W*(T-1), niter);
  } else {
    throw std::runtime_error("mu0==0 not handled yet");
  }

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

    mat mu_t_k;
    if (mu0) {
      mu_t_k.resize(niter, T-1);
      for (int t=1; t<T; t++) {
	for (int i=0; i<niter; i++) {
	  // mu_t_k(i, t-1) = mu_t_vec[k*(T-1)*niter + (t-1)*niter + i];
	  mu_t_k(i, t-1) = mu_t_vec[(t-1)*W*niter + k*niter + i];
	}
      }
    } else {
      throw std::runtime_error("mu0==0 not handled yet");
      // mu_t_k.resize(niter, T-1);
      // for (int t=1; t<T; t++) {
      // 	for (int i=0; i<niter; i++) {
      // 	  mu_t_k(i, t-1) = mu_t_vec[];
      // 	}
      // }
    }

    // for (int t=1; t<3; t++) {
    //   for (int i=0; i<3; i++) {
    // 	std::cout << "k " << k << " t " << t << " i " << i << " " << mu_t_k(i, t-1) << std::endl;
    //   }
    // }

    mat alpha_s(N_knots, niter);
    for (int i=0; i<niter; i++) {
      for (int v=0; v<N_knots; v++) {
        alpha_s(v, i) = alpha_s_vec[v*W*niter + k*niter + i];
      }
    }
    
    cube alpha_t(N_knots, niter, T-1);
    for (int t=1; t<T; t++) {
      for (int i=0; i<niter; i++) {
        for (int v=0; v<N_knots; v++) {
          alpha_t(v, i, t-1) = alpha_t_vec[v*W*(T-1)*niter + k*(T-1)*niter + (t-1)*niter + i];
        }
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
	// for (int j=0; j<N; j++) {
	//   mu_g(j*T, k, i) = mu(i,k) + mu_t(j, k) + Halpha_s(j, k, i);
	// }
	//                mu_g[mu_g_idx,k,i] = mu[i,k] + mu_t_k[i,1] + Halpha_s[,k,i];
      }

      double lambda_inv = 1.0 / lambda(i,k);
      mat Q = exp(-lambda_inv * d_knots);
      mat q = exp(-lambda_inv * d_inter);
      mat Q_inv = inv(Q);
      mat q_Qinv = q * Q_inv;

      vec cs_Csinv_alpha_s = cs_Csinv * alpha_s.col(i);
      // if (i == 3) {
      // 	for (int j=0; j<10; j++) {
      // 	  std::cout << k << " " << j << " " << cs_Csinv_alpha_s(j) << std::endl;
      // 	}
      // }

    
      // if ((i >995) && (i<1000)) {
      // 	for (int t=1; t<3; t++) {
      // 	  // std::cout << i << " " << k << " " << mu(i,k) << std::endl;
      // 	  std::cout << k << " " << i << " " << t-1 << " " << mu_t_k(i,t-1) << std::endl;
      // 	}
      //  }

      for (int t=1; t<T; t++) {

        if (mu0) {
	  vec alpha_t__ = alpha_t(0, i, t-1, size(N_knots, 1, 1));
          Halpha_t(0, k*(T-1) + t-1, i, size(N, 1, 1)) = q_Qinv * alpha_t__;
	  for (int j=0; j<N; j++) {
	    mu_g(j*T+t, k, i) = mu(i,k) + mu_t_k(i, t-1) + cs_Csinv_alpha_s(j) + Halpha_t(j, k*(T-1) + t-1, i);
	  }
        } else {
	  throw std::runtime_error("mu0==0 not handled yet");
	  // vec q_Qinv_alpha_t = q_Qinv * alpha_t(0, i, t-1, size(N_knots, 1, 1));
	  // for (int j=0; j<N; j++) {
	  //   mu_g(j*T+t, k, i) = mu(i,k) + mut_t() + cs_Csinv * alpha_s.col(i) + q_Qinv_alpha_t(j);
	  // }
        }
      }

    }
  }

  return Rcpp::List::create(Rcpp::Named("mu_g") = mu_g,
			    Rcpp::Named("mu") = mu,
			    //			    Rcpp::Named("mu_t") = mu_t,
			    Rcpp::Named("Halpha_s") = Halpha_s,
			    Rcpp::Named("Halpha_t") = Halpha_t);

}
