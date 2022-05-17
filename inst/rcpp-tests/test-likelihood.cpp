#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "../include/intsurv.h"

// check likelihood function of coxph cure rate model for the given new data
// [[Rcpp::export]]
Rcpp::List rt_ll_CoxphCure(const arma::mat& train_surv_x,
                           const arma::mat& train_cure_x,
                           const arma::vec& train_time,
                           const arma::vec& train_event,
                           const arma::mat& test_surv_x,
                           const arma::mat& test_cure_x,
                           const arma::vec& test_time,
                           const arma::vec& test_event,
                           const bool cure_intercept = true,
                           const arma::vec& train_surv_offset = 0,
                           const arma::vec& train_cure_offset = 0,
                           const arma::vec& test_surv_offset = 0,
                           const arma::vec& test_cure_offset = 0,
                           const bool surv_standardize = true,
                           const bool cure_standardize = true,
                           const unsigned int tail_completion = 1,
                           const double tail_tau = -1)
{
    intsurv::Control control0;
    control0.cure(tail_completion, tail_tau);
    intsurv::Control surv_control, cure_control;
    surv_control.set_offset(train_surv_offset)->
        set_standardize(surv_standardize);
    cure_control.logistic(cure_intercept)->
        set_offset(train_cure_offset)->
        set_standardize(cure_standardize);
    // define object
    intsurv::CoxphCure obj {
        train_time, train_event,
        train_surv_x, train_cure_x,
        control0, surv_control, cure_control
    };
    // model-fitting
    obj.fit();
    // set new offsets
    surv_control.set_offset(test_surv_offset);
    cure_control.set_offset(test_cure_offset);
    intsurv::CoxphCure new_obj {
        test_time, test_event,
        test_surv_x, test_cure_x,
        control0, surv_control, cure_control
    };
    double new_ll { obj.obs_log_likelihood(new_obj) };
    return Rcpp::List::create(
        Rcpp::Named("train_negLogL") = obj.neg_ll_,
        Rcpp::Named("test_negLogL") = - new_ll,
        Rcpp::Named("new_surv_x") = new_obj.surv_obj_.get_x(true, false),
        Rcpp::Named("new_cure_x") = new_obj.cure_obj_.get_x(true, false),
        Rcpp::Named("new_time") = new_obj.surv_obj_.time_,
        Rcpp::Named("new_event") = new_obj.surv_obj_.event_,
        Rcpp::Named("new_surv_offset") = intsurv::arma2rvec(
            new_obj.surv_obj_.control_.offset_),
        Rcpp::Named("new_cure_offset") = intsurv::arma2rvec(
            new_obj.cure_obj_.control_.offset_),
        Rcpp::Named("surv_standaridze") = obj.surv_obj_.control_.standardize_,
        Rcpp::Named("cure_standaridze") = obj.cure_obj_.control_.standardize_,
        Rcpp::Named("new_surv_standaridze") =
        new_obj.surv_obj_.control_.standardize_,
        Rcpp::Named("new_cure_standaridze") =
        new_obj.cure_obj_.control_.standardize_,
        Rcpp::Named("ord") = new_obj.surv_obj_.ord_
        );
}

// for Coxph cure rate model with MAR event indicators
// [[Rcpp::export]]
Rcpp::List rt_ll_CoxphCureMar(const arma::mat& train_surv_x,
                              const arma::mat& train_cure_x,
                              const arma::mat& train_mar_x,
                              const arma::vec& train_time,
                              const arma::vec& train_event,
                              const arma::mat& test_surv_x,
                              const arma::mat& test_cure_x,
                              const arma::mat& test_mar_x,
                              const arma::vec& test_time,
                              const arma::vec& test_event,
                              const bool cure_intercept = true,
                              const arma::vec& train_surv_offset = 0,
                              const arma::vec& train_cure_offset = 0,
                              const arma::vec& train_mar_offset = 0,
                              const arma::vec& test_surv_offset = 0,
                              const arma::vec& test_cure_offset = 0,
                              const arma::vec& test_mar_offset = 0,
                              const bool surv_standardize = true,
                              const bool cure_standardize = true,
                              const bool mar_standardize = true,
                              const unsigned int tail_completion = 1,
                              const double tail_tau = -1)
{
    intsurv::Control control0;
    control0.cure(tail_completion, tail_tau);
    intsurv::Control surv_control, cure_control, mar_control;
    surv_control.set_offset(train_surv_offset)->
        set_standardize(surv_standardize);
    cure_control.logistic(cure_intercept)->
        set_offset(train_cure_offset)->
        set_standardize(cure_standardize);
    mar_control.logistic(true)->
        set_offset(train_mar_offset)->
        set_standardize(mar_standardize);
    // define object
    intsurv::CoxphCureMar obj {
        train_time, train_event,
        train_surv_x, train_cure_x, train_mar_x,
        control0, surv_control, cure_control, mar_control
    };
    // model-fitting
    obj.mar_fit();
    obj.fit();
    // update offsets
    surv_control.set_offset(test_surv_offset);
    cure_control.set_offset(test_cure_offset);
    mar_control.set_offset(test_mar_offset);
    intsurv::CoxphCureMar new_obj {
        test_time, test_event,
        test_surv_x, test_cure_x, test_mar_x,
        control0, surv_control, cure_control, mar_control
    };
    double new_ll { obj.obs_log_likelihood(new_obj) };
    return Rcpp::List::create(
        Rcpp::Named("train_negLogL") = obj.neg_ll_,
        Rcpp::Named("test_negLogL") = - new_ll,
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
        Rcpp::Named("ord") = new_obj.surv_obj_.ord_
        );
}
