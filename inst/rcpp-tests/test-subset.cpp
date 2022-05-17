#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "../include/intsurv.h"

// check the subset method for CoxphCure objects
// [[Rcpp::export]]
Rcpp::List rt_subset_CoxphCure(const arma::uvec& index,
                               const arma::mat& surv_x,
                               const arma::mat& cure_x,
                               const arma::vec& time,
                               const arma::vec& event,
                               const bool cure_intercept = true,
                               const arma::vec& surv_offset = 0,
                               const arma::vec& cure_offset = 0,
                               const bool surv_standardize = true,
                               const bool cure_standardize = true)
{
    intsurv::Control control0, surv_control, cure_control;
    surv_control.set_offset(surv_offset)->
        set_standardize(surv_standardize);
    cure_control.logistic(cure_intercept)->
        set_offset(cure_offset)->
        set_standardize(cure_standardize);
    // define object
    intsurv::CoxphCure obj {
        time, event, surv_x, cure_x,
        control0, surv_control, cure_control
    };
    intsurv::CoxphCure new_obj { subset(obj, index) };
    return Rcpp::List::create(
        Rcpp::Named("new_surv_x") = new_obj.surv_obj_.get_x(true, false),
        Rcpp::Named("new_cure_x") = new_obj.cure_obj_.get_x(true, false),
        Rcpp::Named("new_time") = new_obj.surv_obj_.time_,
        Rcpp::Named("new_event") = new_obj.surv_obj_.event_,
        Rcpp::Named("new_surv_offset") = intsurv::arma2rvec(
            new_obj.surv_obj_.control_.offset_),
        Rcpp::Named("new_cure_offset") = intsurv::arma2rvec(
            new_obj.cure_obj_.control_.offset_),
        Rcpp::Named("new_surv_standaridze") =
        new_obj.surv_obj_.control_.standardize_,
        Rcpp::Named("new_cure_standaridze") =
        new_obj.cure_obj_.control_.standardize_,
        Rcpp::Named("rev_ord") = new_obj.surv_obj_.rev_ord_,
        Rcpp::Named("ord") = new_obj.surv_obj_.ord_
        );
}

// for Coxph cure rate model with MAR event indicators
// [[Rcpp::export]]
Rcpp::List rt_subset_CoxphCureMar(const arma::uvec& index,
                                  const arma::mat& surv_x,
                                  const arma::mat& cure_x,
                                  const arma::mat& mar_x,
                                  const arma::vec& time,
                                  const arma::vec& event,
                                  const bool cure_intercept = true,
                                  const bool mar_intercept = true,
                                  const arma::vec& surv_offset = 0,
                                  const arma::vec& cure_offset = 0,
                                  const arma::vec& mar_offset = 0,
                                  const bool surv_standardize = true,
                                  const bool cure_standardize = true,
                                  const bool mar_standardize = true)
{
    intsurv::Control control0, surv_control, cure_control, mar_control;
    surv_control.set_offset(surv_offset)->
        set_standardize(surv_standardize);
    cure_control.logistic(cure_intercept)->
        set_offset(cure_offset)->
        set_standardize(cure_standardize);
    mar_control.logistic(mar_intercept)->
        set_offset(mar_offset)->
        set_standardize(mar_standardize);
    // define object
    intsurv::CoxphCureMar obj {
        time, event, surv_x, cure_x, mar_x,
        control0, surv_control, cure_control, mar_control
    };
    intsurv::CoxphCureMar new_obj { subset(obj, index) };
    return Rcpp::List::create(
        Rcpp::Named("new_surv_x") = new_obj.surv_obj_.get_x(true, false),
        Rcpp::Named("new_cure_x") = new_obj.cure_obj_.get_x(true, false),
        Rcpp::Named("new_mar_x") = new_obj.cure_obj_.get_x(true, false),
        Rcpp::Named("new_time") = new_obj.surv_obj_.time_,
        Rcpp::Named("new_event") = new_obj.surv_obj_.event_,
        Rcpp::Named("new_surv_offset") = intsurv::arma2rvec(
            new_obj.surv_obj_.control_.offset_),
        Rcpp::Named("new_cure_offset") = intsurv::arma2rvec(
            new_obj.cure_obj_.control_.offset_),
        Rcpp::Named("new_mar_offset") = intsurv::arma2rvec(
            new_obj.mar_obj_.control_.offset_),
        Rcpp::Named("new_surv_standaridze") =
        new_obj.surv_obj_.control_.standardize_,
        Rcpp::Named("new_cure_standaridze") =
        new_obj.cure_obj_.control_.standardize_,
        Rcpp::Named("new_mar_standaridze") =
        new_obj.mar_obj_.control_.standardize_,
        Rcpp::Named("rev_ord") = new_obj.surv_obj_.rev_ord_,
        Rcpp::Named("ord") = new_obj.surv_obj_.ord_
        );
}
