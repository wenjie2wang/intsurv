#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "rcpp_coxph.hpp"
#include "l-bfgs.hpp"


// A hopefully fast subroutine fitting Cox model
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_coxph(const arma::vec& time,
                               const arma::vec& event,
                               const arma::mat& z)
{
    Intsurv::RcppCoxph coxph_object {Intsurv::RcppCoxph(time, event, z)};
    arma::vec beta {arma::zeros(z.n_cols)};
    Intsurv::control_lbfgs control {Intsurv::control_lbfgs()};
    Intsurv::control_line_search control_ls {Intsurv::control_line_search()};
    // TODO: may tweak the parameters based on inputs
    control_ls.condition = 2;
    control.m = 8;
    Intsurv::lbfgs(beta, coxph_object, control);
    return Rcpp::NumericVector(beta.begin(), beta.end());
}
