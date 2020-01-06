//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2020  Wenjie Wang <wjwang.stat@gmail.com>
//
// This file is part of the R package intsurv.
//
// The R package intsurv is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package intsurv is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#ifndef SPLINES_H
#define SPLINES_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Intsurv {

    // generate B-spline bases
    inline arma::mat bSpline(const arma::vec& x,
                             const unsigned int& degree,
                             const arma::vec& internal_knots,
                             const arma::vec& boundary_knots)
    {
        const unsigned int order { degree + 1 };
        arma::vec all_knots {
            arma::join_vert(internal_knots, rep_each(boundary_knots, order))
        };
        all_knots = arma::sort(all_knots);

        // degree of freedom
        const unsigned int df { degree + internal_knots.n_elem + 1 };
        arma::mat res { arma::zeros(x.n_elem, df) };

        // piece-wise constant basis
        if (degree == 0) {
            // for each x
            for (size_t i {0}; i < x.n_elem; ++i) {
                // for each basis
                for (size_t j {0}; j < df - 1; ++j) {
                    if (x(i) >= all_knots(j) && x(i) < all_knots(j + 1)) {
                        res(i, j) = 1;
                    }
                }
                // close on right boundary knot
                if (x(i) >= all_knots(df - 1) && x(i) <= all_knots(df)) {
                    res(i, df - 1) = 1;
                }
            }
        } else {
            // use the Cox-de Boor recursive formula
            arma::vec w1 { arma::zeros(x.n_elem) }, w2 { w1 };
            arma::mat b_mat {
                bSpline(x, degree - 1, internal_knots, boundary_knots)
            };
            // for each new basis
            for (size_t j {0}; j < df; ++j) {
                if (isAlmostEqual(all_knots(j + degree), all_knots(j))) {
                    w2 = (all_knots(j + order) - x) /
                        (all_knots(j + order) - all_knots(j + 1));
                    res.col(j) = w2 % b_mat.col(j);
                    continue;
                }
                if (isAlmostEqual(all_knots(j + order), all_knots(j + 1))) {
                    w1 = (x - all_knots(j)) /
                        (all_knots(j + degree) - all_knots(j));
                    res.col(j) = w1 % b_mat.col(j - 1);
                    continue;
                }
                w1 = (x - all_knots(j)) /
                    (all_knots(j + degree) - all_knots(j));
                w2 = (all_knots(j + order) - x) /
                    (all_knots(j + order) - all_knots(j + 1));
                res.col(j) = w1 % b_mat.col(j - 1) + w2 % b_mat.col(j);
            }
        }
        return res;
    }

    // generate derivative of B-splines
    inline arma::mat dbs(const arma::vec& x,
                         const unsigned int& degree,
                         const arma::vec& internal_knots,
                         const arma::vec& boundary_knots,
                         const unsigned int derivs = 1)
    {
        const unsigned int order { degree + 1 };
        const unsigned int df { order + internal_knots.n_elem };
        arma::vec all_knots {
            arma::join_vert(internal_knots, rep_each(boundary_knots, order))
        };
        all_knots = arma::sort(all_knots);

        arma::mat res_mat { arma::zeros(x.n_elem, df) };
        // check for early exit
        if (derivs > degree) {
            return res_mat;
        }
        double double_degree { static_cast<double>(degree) };

        // for first derivative
        arma::mat b_mat {
            bSpline(x, degree - 1, internal_knots, boundary_knots)
        };
        // use a recursive formula
        double w1 { 0 }, w2 { w1 };
        // for higher order derivative
        if (derivs > 1) {
            b_mat = dbs(x, degree - 1, internal_knots,
                        boundary_knots, derivs - 1);
        }
        // for each new basis
        for (size_t j {0}; j < df; ++j) {
            if (isAlmostEqual(all_knots(j + degree), all_knots(j))) {
                w2 = - double_degree /
                    (all_knots(j + order) - all_knots(j + 1));
                res_mat.col(j) = w2 * b_mat.col(j);
                continue;
            }
            if (isAlmostEqual(all_knots(j + order), all_knots(j + 1))) {
                w1 = double_degree /
                    (all_knots(j + degree) - all_knots(j));
                res_mat.col(j) = w1 * b_mat.col(j - 1);
                continue;
            }
            w1 = double_degree /
                (all_knots(j + degree) - all_knots(j));
            w2 = - double_degree /
                (all_knots(j + order) - all_knots(j + 1));
            res_mat.col(j) = w1 * b_mat.col(j - 1) + w2 * b_mat.col(j);
        }
        return res_mat;
    }

    // generate M-spline bases from transformation of B-splines
    inline arma::mat mSpline(const arma::vec& x,
                             const unsigned int& degree,
                             const arma::vec& internal_knots,
                             const arma::vec& boundary_knots)
    {
        const unsigned int order { degree + 1 };
        arma::vec all_knots {
            arma::join_vert(internal_knots, rep_each(boundary_knots, order))
        };
        all_knots = arma::sort(all_knots);
        // B-splines
        arma::mat res_mat {
            bSpline(x, degree, internal_knots, boundary_knots)
        };
        // rescale for each column
        for (size_t j {0}; j < res_mat.n_cols; ++j) {
            res_mat.col(j) *= order / (all_knots(j + order) - all_knots(j));
        }
        return res_mat;
    }

    // generate M-spline bases from transformation of B-splines
    inline arma::mat dms(const arma::vec& x,
                         const unsigned int& degree,
                         const arma::vec& internal_knots,
                         const arma::vec& boundary_knots,
                         const unsigned int derivs = 1)
    {
        const unsigned int order { degree + 1 };
        arma::vec all_knots {
            arma::join_vert(internal_knots, rep_each(boundary_knots, order))
        };
        all_knots = arma::sort(all_knots);
        // derivatives of B-splines
        arma::mat res_mat {
            dbs(x, degree, internal_knots, boundary_knots, derivs)
        };
        // rescale for each column
        for (size_t j {0}; j < res_mat.n_cols; ++j) {
            res_mat.col(j) *= order / (all_knots(j + order) - all_knots(j));
        }
        return res_mat;
    }

    // generate I-splines from B-splines
    inline arma::mat iSpline(const arma::vec& x,
                             const unsigned int& degree,
                             const arma::vec& internal_knots,
                             const arma::vec& boundary_knots)
    {
        const unsigned int n_knots { internal_knots.n_elem };
        const unsigned int order { degree + 1 };
        const unsigned int df { order + n_knots };

        // determine j from x
        arma::vec j_vec { rep_double(order, x.n_elem) };
        if (internal_knots.n_elem > 0) {
            j_vec = step_fun(x, internal_knots,
                             arma::regspace<arma::vec>(order, df));
        }
        // generate B-splines
        arma::mat bs_mat {
            bSpline(x, order, internal_knots, boundary_knots)
                };
        arma::mat res_mat { arma::zeros(x.n_elem, df) };
        // for each x
        for (size_t k {0}; k < x.n_elem; ++k) {
            unsigned int j { static_cast<unsigned int>(j_vec(k)) };
            arma::rowvec bs_row { bs_mat.row(k) };
            // for each basis
            for (size_t i {0}; i < df; ++i) {
                if (i + 1 > j) {
                    res_mat(k, i) = 0;
                } else if (i + 1 < j - degree) {
                    res_mat(k, i) = 1;
                } else {
                    res_mat(k, i) = arma::sum(bs_row.subvec(i + 1, j));
                }
            }
        }
        return res_mat;
    }

    // get the placement of knots based on quantile
    inline arma::vec get_boundary_knots(const arma::vec& x)
    {
        arma::vec probs { {0, 1} };
        return arma_quantile(x, probs);
    }
    inline arma::vec get_internal_knots(const arma::vec& x,
                                        const unsigned int& num_knots)
    {
        double prob { 1.0 / (static_cast<double>(num_knots) + 1.0) };
        arma::vec probs { prob * arma::regspace(1, num_knots) };
        return arma_quantile(x, probs);
    }

    // non-negative least square
    inline arma::vec nnls(const arma::vec& y, const arma::mat& x,
                          const double eps = 1e-4)
    {
        // reference: Lawson and Hanson (1974)
        const unsigned int p { x.n_cols };
        arma::uvec setP { arma::zeros<arma::uvec>(p) };
        arma::uvec setR { arma::ones<arma::uvec>(p) };
        arma::vec beta { arma::zeros(p) };
        arma::vec w { mat2vec(x.t() * y) };

        // main loop
        while (arma::sum(setR) > 0 && w.max() > eps) {
            arma::uvec idxR { arma_which(setR) };
            size_t j { idxR(w.rows(idxR).index_max()) };
            setP(j) = 1;
            setR(j) = 0;
            arma::uvec idxP { arma_which(setP) };
            arma::mat xP { x.cols(idxP) };
            arma::vec s { arma::zeros(p) };
            arma::vec sP {
                mat2vec(arma::inv_sympd(xP.t() * xP) * (xP.t() * y))
            };
            s.rows(idxP) = sP;
            // inner loop
            while (sP.min() <= 0) {
                arma::uvec np_idx { arma::find(sP <= 0) };
                arma::vec betaP { beta.elem(idxP) };
                double alpha {
                    arma::min(betaP.elem(np_idx) /
                              (betaP.elem(np_idx) - sP.elem(np_idx)))
                };
                beta += alpha * (s - beta);
                // move j in setP where beta(j) = 0 to setR
                for (size_t i {0}; i < p; ++i) {
                    if (setP(i) && isAlmostEqual(beta(i), 0)) {
                        setP(i) = 0;
                        setR(i) = 1;
                    }
                }
                idxP = arma_which(setP);
                xP = x.cols(idxP);
                sP = mat2vec(arma::inv_sympd(xP.t() * xP) * (xP.t() * y));
                s = arma::zeros(p);
                s.rows(idxP) = sP;
            }
            beta = s;
            w = mat2vec(x.t() * (y - x * beta));
        }
        return beta;
    }


}

#endif
