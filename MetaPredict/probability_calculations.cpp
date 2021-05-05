//LogL.betabinomial
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logBetabinomial(NumericVector params, NumericVector yl, double n, NumericVector ml) {
  double alpha = params[0];
  double beta = params[1];
  return sum(-lgamma(alpha + beta) - lgamma(alpha + yl) - lgamma(beta + (n * ml) - yl) + lgamma(alpha) + lgamma(beta) + lgamma(alpha + beta + (n * ml)));
}



//get_pj.n_betabinomial.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double get_pjn_betabinomial(double pxj_n, double  pxj_nbar, double alpha, double beta,
                            double yj, double  n, double mj, double xj) {
  double k = n;
  double Aj = alpha + yj;
  double Bj = beta + (n * mj) - yj;
  return exp(lgamma(Aj + Bj) + lgamma(Aj + k) + lgamma(Bj + n - k) - lgamma(Aj) - lgamma(Bj) - lgamma(Aj + Bj + n));
}



//get_pj.k_betabinomial.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_pjk_betabinomial(double pxj_n, double  pxj_nbar, double alpha, double beta,
                            double yj, double  n, double mj, double xj) {
  NumericVector k = seq(xj, (n - 1));
  double Aj = alpha + yj;
  double Bj = beta + (n * mj) - yj;
  return exp(lgamma(Aj + Bj) + lgamma(Aj + k) + lgamma(Bj + n - k) - lgamma(Aj) - lgamma(Bj) - lgamma(Aj + Bj + n));
}







