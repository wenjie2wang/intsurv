#ifndef LOGISTIC_REG_H
#define LOGISTIC_REG_H

#include <RcppArmadillo.h>
#include "utils.hpp"

namespace Intsurv {

    // define class for inputs and outputs
    class LogisticReg
    {
    private:
        arma::mat x;
        arma::vec y;

    public:
        // constructors
        LogisticReg(const arma::mat& x_, const arma::vec& y_) :
            x(x_), y(y_) {}

        // function members
        inline double objective(const arma::vec& beta) const;
        inline arma::vec gradient(const arma::vec& beta) const;
        inline double objective(const arma::vec& beta, arma::vec& grad) const;

    };

    // define objective function
    inline double LogisticReg::objective(const arma::vec& beta) const
    {
        arma::vec x_beta {x * beta};
        arma::vec exp_x_beta {arma::exp(x_beta)};
        arma::vec y_x_beta {y % x_beta};
        double res {
            arma::as_scalar(arma::sum(arma::log(1 + exp_x_beta) - y_x_beta))
        };
        return res;
    }

    // define gradient function
    inline arma::vec LogisticReg::gradient(const arma::vec& beta) const
    {
        arma::vec x_beta {x * beta};
        arma::vec exp_x_beta {arma::exp(x_beta)};
        arma::vec res
            { x.t() * (exp_x_beta / (1 + exp_x_beta) - y)};
        return res;
    }

    // define objective function and overwrites graidient
    inline double LogisticReg::objective(const arma::vec& beta,
                                         arma::vec& grad) const
    {
        arma::vec x_beta {x * beta};
        arma::vec exp_x_beta {arma::exp(x_beta)};
        grad = x.t() * (exp_x_beta / (1 + exp_x_beta) - y);
        arma::vec y_x_beta {y % x_beta};
        double res {
            arma::as_scalar(arma::sum(arma::log(1 + exp_x_beta) - y_x_beta))
        };
        return res;
    }

}


#endif
