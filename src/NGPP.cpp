#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;


// expectations of log(cosh(Z)) and -exp(-Z^2/2) for Z ~ N(0,1)
static const double Elogcosh = 0.374567207491438;
static const double Egauss = -0.707106781186548;

// numerically stable log(cosh(x)) = |x| + log(1 + exp(-2|x|)) - log(2),
// log(cosh(x)) itself overflows for |x| > 710
static arma::mat logcosh(const arma::mat& x) {
  arma::mat ax = abs(x);
  return ax + log(1.0 + exp(-2.0*ax)) - std::log(2.0);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat symmetricPower_C(arma::mat x, double r) {
  
  arma::mat eig_vec;
  arma::vec eig_val;
  
  eig_sym(eig_val, eig_vec, x);
  
  arma::mat pow_val = diagmat(pow(sort(eig_val, "descend"), r));
  
  return fliplr(eig_vec)*pow_val*fliplr(eig_vec).t();
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec computeObj_C(arma::mat x, arma::vec nl, arma::vec alpha) {
  
  int nl_num = nl.n_elem;
  
  arma::vec obj = zeros<arma::vec>(x.n_cols);

  
  for(int i = 0; i < nl_num; i++){
    if(nl(i) == 1){
        obj = obj + alpha(i)*pow(mean(pow(x, 3), 0).t(), 2);
    }
    else if(nl(i) == 2){
      obj = obj + alpha(i)*pow(mean(pow(x, 4), 0).t() - 3, 2);
    }
    else if(nl(i) == 3){
      obj = obj + alpha(i)*pow(mean(logcosh(x), 0).t() - Elogcosh, 2);
    }
    else if(nl(i) == 4){
      obj = obj + alpha(i)*pow(mean(-1*exp(-0.5*pow(x, 2)), 0).t() - Egauss, 2);
    }
  }
  
  return obj;
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec computeTVec_C(arma::vec u, arma::mat x, arma::vec nl, arma::vec alpha) {
  
  int nl_num = nl.n_elem;
  
  arma::vec comp = x*u;
  
  arma::vec T = zeros<arma::vec>(x.n_cols);
  
  
  for(int i = 0; i < nl_num; i++){
    if(nl(i) == 1){
      T = T + alpha(i)*as_scalar(mean(pow(comp, 3), 0))*(3*mean(x.each_col() % pow(comp, 2), 0).t());
    }
    if(nl(i) == 2){
      T = T + alpha(i)*as_scalar(mean(pow(comp, 4), 0) - 3)*(4*mean(x.each_col() % pow(comp, 3), 0).t() - 12*u);
    }
    if(nl(i) == 3){
      T = T + alpha(i)*(as_scalar(mean(logcosh(comp))) - Elogcosh)*(mean(x.each_col() % tanh(comp), 0).t() - as_scalar(1 - mean(pow(tanh(comp), 2)))*u);
    }
    if(nl(i) == 4){
      T = T + alpha(i)*(as_scalar(mean(-1*exp(-0.5*pow(comp, 2)))) - Egauss)*(mean(x.each_col() % (comp % exp(-0.5*pow(comp, 2))), 0).t() - as_scalar(mean((1 - pow(comp, 2)) % exp(-0.5*pow(comp, 2))))*u);
    }
  }
  
  return T;
}



