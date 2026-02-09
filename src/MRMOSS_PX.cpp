#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <set>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <Rcpp.h>
#include <ctime>


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double loglikelihood_cpp(arma::vec gamma_hat, arma::mat Gamma_hat, arma::mat R, int n1, int n2, arma::vec theta){
  int p = Gamma_hat.n_rows;
  int m = Gamma_hat.n_cols;
  
  // Extract parameters from theta vector
  arma::vec beta = theta.subvec(0, m-1);
  double sigma2 = theta(m);
  arma::vec tau2 = theta.subvec(m+1, 2*m);
  double sigmax2 = theta(2*m + 1);
  arma::vec sigma = theta.subvec(2*m + 2, 3*m + 1);
  
  // Construct matrix A
  arma::mat A = arma::eye<arma::mat>(m+1, m+1);
  A.col(0).rows(1, m) = beta;
  
  // Construct matrix Sigma
  arma::mat Sigma = arma::zeros<arma::mat>(m+1, m+1);
  Sigma(0,0) = sigma2;
  Sigma.submat(1,1,m,m) = arma::diagmat(tau2);
  
  // Construct matrix D and Lambda
  arma::mat D = arma::diagmat(sigma);
  arma::mat Lambda = D * R * D;
  
  // Construct matrix Omega
  arma::mat Omega = arma::zeros<arma::mat>(m+1, m+1);
  Omega(0,0) = sigmax2 / n1;
  Omega.submat(1,1,m,m) = Lambda / n2;
  
  // Combine gamma_hat and Gamma_hat into matrix g
  arma::mat g = arma::join_horiz(gamma_hat, Gamma_hat);
  
  // Compute Omega0
  arma::mat Omega0 = A * Sigma * A.t() + Omega;
  
  // Compute log-likelihood components
  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, Omega0);
  
  double term1 = -p * (m+1) * std::log(2 * M_PI) / 2.0;
  double term2 = -p * log_det_val / 2.0;
  double term3 = -arma::trace(arma::inv(Omega0) * g.t() * g) / 2.0;
  
  return term1 + term2 + term3;
}

// [[Rcpp::export]]
List MRMOSS_PX_cpp(arma::vec gamma_hat, arma::mat Gamma_hat, arma::mat R, 
                  int n1, int n2, arma::vec theta0, 
                  arma::uvec test, int maxiter, double rd) {
    // ==== shrinkage correlation matrix ====


    int p = Gamma_hat.n_rows;
    int m = Gamma_hat.n_cols;
    arma::mat g = join_horiz(gamma_hat, Gamma_hat);
    arma::mat R_inv = inv(R);
    
    // Construct list of Vi
    std::vector<arma::mat> V(m);
    for(int i = 0; i < m; i++) {
        arma::mat Vi = zeros<mat>(m+1, m+1);
        Vi(i+1, 0) = 1; 
        V[i] = Vi;
    }
    
    // initialize parameters
    arma::vec theta_new; 
    arma::vec beta_new = theta0.subvec(0, m-1);
    double sigma2_new = theta0(m);
    arma::vec tau2_new = theta0.subvec(m+1, 2*m);
    double sigmax2_new = theta0(2*m + 1);
    arma::vec sigma_new = theta0.subvec(2*m + 2, 3*m + 1);
    
    double loglikeli = 0;
    double dif = 1;
    int iter = 0;
    arma::vec loglikeli_vec(maxiter + 2); 
    loglikeli_vec(0) = 0;
    
    while((iter <= maxiter) && (dif >= 1e-8)) {
        arma::vec beta_old = beta_new;
        double sigma2_old = sigma2_new;
        arma::vec tau2_old = tau2_new;
        double sigmax2_old = sigmax2_new;
        arma::vec sigma_old = sigma_new;
        
        arma::mat Sigma = zeros<mat>(m+1, m+1);
        Sigma(0,0) = sigma2_old;
        Sigma.submat(1, 1, m, m) = diagmat(tau2_old);
        
        arma::mat Sigma_inv = zeros<mat>(m+1, m+1);
        Sigma_inv(0,0) = 1.0 / sigma2_old;
        Sigma_inv.submat(1, 1, m, m) = diagmat(1.0 / tau2_old);
        
        arma::mat A = eye<mat>(m+1, m+1);
        A.col(0).subvec(1, m) = beta_old;
        
        arma::mat D = diagmat(sigma_old);
        arma::mat Lambda = D * R * D;
        arma::mat Lambda_inv = diagmat(1.0 / sigma_old) * R_inv * diagmat(1.0 / sigma_old);
        
        arma::mat Omega = zeros<mat>(m+1, m+1);
        Omega(0,0) = sigmax2_old / n1;
        Omega.submat(1, 1, m, m) = Lambda / n2;
        
        arma::mat Omega_inv = zeros<mat>(m+1, m+1);
        Omega_inv(0,0) = n1 / sigmax2_old;
        Omega_inv.submat(1, 1, m, m) = n2 * Lambda_inv;
        
        // E step
        arma::mat Sigma_tilde_inv = Sigma_inv + A.t() * Omega_inv * A;
        arma::mat Sigma_tilde = inv(Sigma_tilde_inv);
        arma::mat mu0 = Sigma_tilde * A.t() * Omega_inv;
        arma::mat mu = g * mu0.t();
        
        // M step 
        // update lambda
        double lambda = sum(gamma_hat % mu.col(0)) / 
                       (p * Sigma_tilde(0,0) + sum(square(mu.col(0))));
        
        arma::mat A_lambda = A;
        A_lambda(0, 0) = lambda;
        arma::mat G_lambda = A_lambda * Sigma_tilde * A_lambda.t();
        
        // update beta, tau2, sigma
        for(int i = 0; i < m; i++) {
            double beta_num = -p * trace(V[i].t() * Omega_inv * (A_lambda - beta_old(i) * V[i]) * Sigma_tilde);
            double beta_den = p * trace(V[i].t() * Omega_inv * V[i] * Sigma_tilde);
            
            double Ui = 0;
            double Vi = 0;
            
            for(int j = 0; j < p; j++) {
                arma::rowvec g_j = g.row(j);
                arma::vec muj = mu0 * g_j.t();
                
                arma::rowvec temp = muj.t() * V[i].t() * Omega_inv;
                beta_num += as_scalar(temp * (g_j.t() - (A_lambda - beta_old(i) * V[i]) * muj));
                beta_den += as_scalar(temp * V[i] * muj);
                
                arma::vec temp2 = g_j.t() - A_lambda * muj;
                arma::mat Gj_lambda = G_lambda + temp2 * temp2.t();
                
                arma::rowvec Gj_row = Gj_lambda.submat(1, i+1, m, i+1).t();
                Ui += n2 * (as_scalar(Gj_row * (R_inv.row(i).t() / sigma_old)) - 
                          Gj_lambda(i+1, i+1) * R_inv(i,i) / sigma_old(i));
                Vi += n2 * Gj_lambda(i+1, i+1) * R_inv(i,i);
            }
            
            beta_new(i) = beta_num / beta_den;
            tau2_new(i) = Sigma_tilde(i+1, i+1) + sum(square(mu.col(i+1))) / p;
            sigma_new(i) = 2 * Vi / (sqrt(Ui*Ui + 4 * p * Vi) - Ui);
            sigma_new(i) *= rd;
        }
        
        // update sigma2, sigmax2
        sigma2_new = Sigma_tilde(0,0) + sum(square(mu.col(0))) / p;
        sigmax2_new = n1 * lambda * Sigma_tilde(0,0) + n1 * sum(square(gamma_hat - lambda * mu.col(0))) / p;
        
        // Reduction step
        beta_new = beta_new / lambda;
        sigma2_new = lambda * lambda * sigma2_new;
        lambda = 1;
        
        theta_new = join_vert(
          join_vert(
            join_vert(
              join_vert(beta_new, vec({sigma2_new})),
              tau2_new
            ),
            vec({sigmax2_new})
          ),
          sigma_new
        );
        
        loglikeli = loglikelihood_cpp(gamma_hat, Gamma_hat, R, n1, n2, theta_new);
        loglikeli_vec(iter + 1) = loglikeli;
        dif = std::abs(loglikeli_vec(iter + 1) - loglikeli_vec(iter));
        iter++;
    }
    
    // Estimate under null hypothesis Hi0
    arma::mat loglikeli_null = arma::zeros<arma::mat>(test.n_elem, maxiter + 2);
    arma::mat theta_null = arma::zeros<arma::mat>(test.n_elem, theta0.n_elem);
    arma::vec LRT = arma::zeros<arma::vec>(test.n_elem);
    arma::vec pvalue = arma::zeros<arma::vec>(test.n_elem);
    
    for (arma::uword t_idx = 0; t_idx < test.n_elem; ++t_idx) {
      arma::uword t = test(t_idx) - 1;  
      double dif = 1.0;
      int iter0 = 0;
      
      while ((iter0 <= maxiter) && (dif >= 1e-8)) {
        arma::vec beta_old = beta_new;
        beta_old(t) = 0.0;
        double sigma2_old = sigma2_new;
        arma::vec tau2_old = tau2_new;
        double sigmax2_old = sigmax2_new;
        arma::vec sigma_old = sigma_new;
        
        arma::mat Sigma = arma::zeros<arma::mat>(m+1, m+1);
        Sigma(0,0) = sigma2_old;
        Sigma.submat(1,1,m,m) = arma::diagmat(tau2_old);
        
        arma::mat Sigma_inv = arma::zeros<arma::mat>(m+1, m+1);
        Sigma_inv(0,0) = 1.0 / sigma2_old;
        Sigma_inv.submat(1,1,m,m) = arma::diagmat(1.0 / tau2_old);
        
        arma::mat A = arma::eye<arma::mat>(m+1, m+1);
        A.col(0).subvec(1,m) = beta_old;
        
        arma::mat D = arma::diagmat(sigma_old);
        arma::mat Lambda = D * R * D;
        arma::mat Lambda_inv = arma::diagmat(1.0 / sigma_old) * R_inv * arma::diagmat(1.0 / sigma_old);
        
        arma::mat Omega = arma::zeros<arma::mat>(m+1, m+1);
        Omega(0,0) = sigmax2_old / n1;
        Omega.submat(1,1,m,m) = Lambda / n2;
        
        arma::mat Omega_inv = arma::zeros<arma::mat>(m+1, m+1);
        Omega_inv(0,0) = n1 / sigmax2_old;
        Omega_inv.submat(1,1,m,m) = n2 * Lambda_inv;
        
        // E-step
        arma::mat Sigma_tilde_inv = Sigma_inv + A.t() * Omega_inv * A;
        arma::mat Sigma_tilde = arma::inv(Sigma_tilde_inv);
        arma::mat mu0 = Sigma_tilde * A.t() * Omega_inv;
        arma::mat mu = g * mu0.t();
        arma::mat G = A * Sigma_tilde * A.t();
        
        // M-step
        for (arma::uword i = 0; i < m; ++i) {
          double beta_num = -p * arma::trace(V[i].t() * Omega_inv * (A - beta_old(i) * V[i]) * Sigma_tilde);
          double beta_den = p * arma::trace(V[i].t() * Omega_inv * V[i] * Sigma_tilde);
          
          double Ui = 0.0;
          double Vi = 0.0;
          
          for (arma::uword j = 0; j < p; ++j) {
            arma::rowvec g_j = g.row(j);
            arma::vec muj = mu0 * g_j.t();
            arma::rowvec temp = muj.t() * V[i].t() * Omega_inv;
            
            if (i == t) {
              beta_num = 0.0;
              beta_den = 1.0;
            } else {
              beta_num += arma::as_scalar(temp * (g_j.t() - (A - beta_old(i) * V[i]) * muj));
              beta_den += arma::as_scalar(temp * V[i] * muj);
            }
            
            arma::vec temp2 = g_j.t() - A * muj;
            arma::mat Gj = G + temp2 * temp2.t();
            
            arma::rowvec Gj_row = Gj.submat(1, i+1, m, i+1).t();
            Ui += n2 * (arma::as_scalar(Gj_row * (R_inv.row(i).t() / sigma_old)) - 
              Gj(i+1, i+1) * R_inv(i,i) / sigma_old(i));
            Vi += n2 * Gj(i+1, i+1) * R_inv(i,i);
          }
          
          beta_new(i) = beta_num / beta_den;
          tau2_new(i) = Sigma_tilde(i+1, i+1) + arma::sum(arma::square(mu.col(i+1))) / p;
          sigma_new(i) = 2.0 * Vi / (std::sqrt(Ui*Ui + 4.0 * p * Vi) - Ui);
          sigma_new(i) *= rd;
        }
        
        sigma2_new = Sigma_tilde(0,0) + arma::sum(arma::square(mu.col(0))) / p;
        sigmax2_new = n1 * Sigma_tilde(0,0) + n1 * arma::sum(arma::square(gamma_hat - mu.col(0))) / p;
        
        theta_null.row(t_idx) = arma::join_vert(
          arma::join_vert(
            arma::join_vert(
              arma::join_vert(beta_new, arma::vec({sigma2_new})),
              tau2_new
            ),
            arma::vec({sigmax2_new})
          ),
          sigma_new
        ).t();
        
        loglikeli_null(t_idx, iter0 + 1) = loglikelihood_cpp(gamma_hat, Gamma_hat, R, n1, n2, theta_null.row(t_idx).t());
        
        ++iter0;
        dif = std::abs(loglikeli_null(t_idx, iter0) - loglikeli_null(t_idx, iter0 - 1));
      }
      
      LRT(t_idx) = 2.0 * (loglikeli_vec(iter) - loglikeli_null(t_idx, iter0));
      pvalue(t_idx) = R::pchisq(LRT(t_idx), 1.0, false, false);
    }
    
    // Estimate under overall null hypothesis H0
    arma::mat loglikeli_overall = arma::zeros<arma::mat>(1, maxiter + 2);
    arma::mat theta_overall = arma::zeros<arma::mat>(1, theta0.n_elem);
    double LRT_overall = 0;
    double pvalue_overall = 0;
    
    double dif_overall = 1.0;
    int iter_overall = 0;
    
    while ((iter_overall <= maxiter) && (dif_overall >= 1e-8)) {
      arma::vec beta_old = arma::zeros<arma::vec>(m);
      arma::vec beta_new = beta_old;
      double sigma2_old = sigma2_new;
      arma::vec tau2_old = tau2_new;
      double sigmax2_old = sigmax2_new;
      arma::vec sigma_old = sigma_new;
      
      arma::mat Sigma = arma::zeros<arma::mat>(m+1, m+1);
      Sigma(0,0) = sigma2_old;
      Sigma.submat(1,1,m,m) = arma::diagmat(tau2_old);
      
      arma::mat Sigma_inv = arma::zeros<arma::mat>(m+1, m+1);
      Sigma_inv(0,0) = 1.0 / sigma2_old;
      Sigma_inv.submat(1,1,m,m) = arma::diagmat(1.0 / tau2_old);
      
      arma::mat A = arma::eye<arma::mat>(m+1, m+1);
      A.col(0).subvec(1,m) = beta_old;
      
      arma::mat D = arma::diagmat(sigma_old);
      arma::mat Lambda = D * R * D;
      arma::mat Lambda_inv = arma::diagmat(1.0 / sigma_old) * R_inv * arma::diagmat(1.0 / sigma_old);
      
      arma::mat Omega = arma::zeros<arma::mat>(m+1, m+1);
      Omega(0,0) = sigmax2_old / n1;
      Omega.submat(1,1,m,m) = Lambda / n2;
      
      arma::mat Omega_inv = arma::zeros<arma::mat>(m+1, m+1);
      Omega_inv(0,0) = n1 / sigmax2_old;
      Omega_inv.submat(1,1,m,m) = n2 * Lambda_inv;
      
      // E-step
      arma::mat Sigma_tilde_inv = Sigma_inv + A.t() * Omega_inv * A;
      arma::mat Sigma_tilde = arma::inv(Sigma_tilde_inv);
      arma::mat mu0 = Sigma_tilde * A.t() * Omega_inv;
      arma::mat mu = g * mu0.t();
      arma::mat G = A * Sigma_tilde * A.t();
      
      // M-step
      for (arma::uword i = 0; i < m; ++i) {
        double Ui = 0.0;
        double Vi = 0.0;
        
        for (arma::uword j = 0; j < p; ++j) {
          arma::rowvec g_j = g.row(j);
          arma::vec muj = mu0 * g_j.t();
          arma::rowvec temp = muj.t() * V[i].t() * Omega_inv;
          
          arma::vec temp2 = g_j.t() - A * muj;
          arma::mat Gj = G + temp2 * temp2.t();
          
          arma::rowvec Gj_row = Gj.submat(1, i+1, m, i+1).t();
          Ui += n2 * (arma::as_scalar(Gj_row * (R_inv.row(i).t() / sigma_old)) - 
            Gj(i+1, i+1) * R_inv(i,i) / sigma_old(i));
          Vi += n2 * Gj(i+1, i+1) * R_inv(i,i);
        }
        tau2_new(i) = Sigma_tilde(i+1, i+1) + arma::sum(arma::square(mu.col(i+1))) / p;
        sigma_new(i) = 2.0 * Vi / (std::sqrt(Ui*Ui + 4.0 * p * Vi) - Ui);
        sigma_new(i) *= rd;
      }
      
      sigma2_new = Sigma_tilde(0,0) + arma::sum(arma::square(mu.col(0))) / p;
      sigmax2_new = n1 * Sigma_tilde(0,0) + n1 * arma::sum(arma::square(gamma_hat - mu.col(0))) / p;
      
      theta_overall = arma::join_vert(
        arma::join_vert(
          arma::join_vert(
            arma::join_vert(beta_new, arma::vec({sigma2_new})),
            tau2_new
          ),
          arma::vec({sigmax2_new})
        ),
        sigma_new
      );
      
      loglikeli_overall(iter_overall + 1) = loglikelihood_cpp(gamma_hat, Gamma_hat, R, n1, n2, theta_overall);
      
      ++iter_overall;
      dif_overall = std::abs(loglikeli_overall(iter_overall) - loglikeli_overall(iter_overall - 1));
    }
    
    LRT_overall = 2.0 * (loglikeli_vec(iter) - loglikeli_overall(iter_overall));
    pvalue_overall = R::pchisq(LRT_overall, m, false, false);
    
    return Rcpp::List::create(
      Rcpp::Named("beta") = theta_new.subvec(0, m-1).t(),
      Rcpp::Named("theta") = theta_new.t(),
      Rcpp::Named("theta_null") = theta_null,
      Rcpp::Named("loglikeli") = loglikeli_vec.subvec(1, iter),
      Rcpp::Named("loglikeli_null") = loglikeli_null.cols(1, loglikeli_null.n_cols - 1),
      Rcpp::Named("iteration") = iter - 1,
      Rcpp::Named("LRT") = LRT,
      Rcpp::Named("pvalue") = pvalue,
      Rcpp::Named("theta_overall") = theta_overall,
      Rcpp::Named("pvalue_overall") = pvalue_overall
    );
}

