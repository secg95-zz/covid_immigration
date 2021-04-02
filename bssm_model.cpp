// Model definition for bssm package

// #include <armadillo>
// #include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// Function for the prior mean of alpha_1
// [[Rcpp::export]]
arma::vec a1_fn(const arma::vec& theta, const arma::vec& known_params) {
  arma::vec a1(3);
  a1(0) = known_params(0);
  a1(1) = known_params(1);
  a1(2) = known_params(2);
  return a1;
}

// Function for the prior covariance matrix of alpha_1
// [[Rcpp::export]]
arma::mat P1_fn(const arma::vec& theta, const arma::vec& known_params) {
  arma::mat P1(3, 3);
  P1(0, 0) = known_params(3);
  P1(1, 1) = known_params(4);
  P1(2, 2) = known_params(5);
  P1(0, 1) = known_params(6);
  P1(1, 0) = known_params(6);
  P1(0, 2) = known_params(7);
  P1(2, 0) = known_params(7);
  P1(1, 2) = known_params(8);
  P1(2, 1) = known_params(8);
  return P1;
}

arma::mat T_prime = {
  {1, 0, 1},
  {0, 1, 0},
  {0, 0, 1}
};

// [[Rcpp::export]]
arma::vec T_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params) {
  return T_prime * alpha;
}

// [[Rcpp::export]]
arma::mat T_gn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params) {
  return T_prime;
}

// [[Rcpp::export]]
arma::mat R_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat R(3,3); 
  R.eye();
  return R;
}

// [[Rcpp::export]]
arma::mat H_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat H(3,3); 
  H.zeros();
  return H;
}

// [[Rcpp::export]]
arma::vec Z_fn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params) {
  // theta := infectivity
  arma::vec ans(1);
  ans(0) = exp(alpha(0)) * (theta(t) + exp(alpha(1)));
  return ans;
}

// [[Rcpp::export]]
arma::mat Z_gn(const unsigned int t, const arma::vec& alpha, const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params) {
  arma::mat ans(1, 3);
  ans(0, 0) = exp(alpha(0)) * (theta(t) + exp(alpha(1)));
  ans(0, 1) = exp(alpha(0)) * exp(alpha(1));
  ans(0, 2) = 0;
  return ans;
}

// [[Rcpp::export]]
double log_prior_pdf(const arma::vec& theta) {
  return 1;
}

// [[Rcpp::export]]
Rcpp::List create_xptrs() {
  // typedef for a pointer of nonlinear function of model equation returning vec (T, Z)
  typedef arma::vec (*nvec_fnPtr)(const unsigned int t, const arma::vec& alpha, 
                     const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params);
  // typedef for a pointer of nonlinear function returning mat (Tg, Zg, H, R)
  typedef arma::mat (*nmat_fnPtr)(const unsigned int t, const arma::vec& alpha, 
                     const arma::vec& theta, const arma::vec& known_params, const arma::mat& known_tv_params);
  // typedef for a pointer returning a1
  typedef arma::vec (*a1_fnPtr)(const arma::vec& theta, const arma::vec& known_params);
  // typedef for a pointer returning P1
  typedef arma::mat (*P1_fnPtr)(const arma::vec& theta, const arma::vec& known_params);
  // typedef for a pointer of log-prior function
  typedef double (*prior_fnPtr)(const arma::vec& theta);
  
  return Rcpp::List::create(
    Rcpp::Named("a1_fn") = Rcpp::XPtr<a1_fnPtr>(new a1_fnPtr(&a1_fn)),
    Rcpp::Named("P1_fn") = Rcpp::XPtr<P1_fnPtr>(new P1_fnPtr(&P1_fn)),
    Rcpp::Named("Z_fn") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&Z_fn)),
    Rcpp::Named("H_fn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&H_fn)),
    Rcpp::Named("T_fn") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&T_fn)),
    Rcpp::Named("R_fn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&R_fn)),
    Rcpp::Named("Z_gn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&Z_gn)),
    Rcpp::Named("T_gn") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&T_gn)),
    Rcpp::Named("log_prior_pdf") = 
      Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf)));
}

/*
int main() {
  arma::vec alpha = {1, 1, 1};
  arma::vec theta = {1, 2, 3, 4, 5, 6, 7, 8};
  // std::cout << T_fn(0, alpha, arma::vec(), arma::vec(), arma::mat()) << std::endl;
  std::cout << Z_fn(3, alpha, theta, arma::vec(), arma::mat()) << std::endl;
}
*/