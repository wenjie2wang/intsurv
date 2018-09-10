#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "logistic_reg.hpp"
#include "l-bfgs.hpp"


// A hopefully fast subroutine fitting logistic model
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_logistic(const arma::mat& x, const arma::vec y)
{
    Intsurv::LogisticReg object {Intsurv::LogisticReg(x, y)};
    arma::vec beta {arma::zeros(x.n_cols)};
    Intsurv::control_line_search control_ls {
        Intsurv::control_line_search()
    };
    Intsurv::control_lbfgs control {Intsurv::control_lbfgs()};
    // TODO: may tweak the parameters based on inputs
    control.line_search = control_ls;
    control.epsilon = 1e-3;
    control.m = 6;
    Intsurv::lbfgs(beta, object, control);
    return Rcpp::NumericVector(beta.begin(), beta.end());
}
