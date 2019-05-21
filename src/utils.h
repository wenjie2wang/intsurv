//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2019  Wenjie Wang <wjwang.stat@gmail.com>
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

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>            // std::max, std::set_union, etc.
#include <cmath>                // std::pow and std::sqrt, etc.
#include <limits>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

// #include <armadillo>
#include <RcppArmadillo.h>

namespace Intsurv {

    // compare double-precision numbers for almost equality
    inline bool isAlmostEqual(double A, double B)
    {
        double MaxRelDiff {std::numeric_limits<double>::epsilon()};
        // compute the difference.
        double diff = std::abs(A - B);
        A = std::abs(A);
        B = std::abs(B);
        // Find the largest
        double largest = (B > A) ? B : A;
        if (diff <= largest * MaxRelDiff) {
            return true;
        } else {
            return false;
        }
    }

    // function checking if there exists any duplicates
    inline bool any_duplicated(const arma::vec& x)
    {
        std::unordered_set<double> seen;
        bool res {false};
        for (size_t i {0}; i < x.n_rows; ++i) {
            res = ! seen.insert(x(i)).second;
            if (res) break;
        }
        return res;
    }
    // returns indices of duplicated elements
    inline arma::uvec duplicated(const arma::vec& x, bool fromLast = false)
    {
        std::unordered_set<double> seen;
        std::vector<unsigned int> res;
        if (fromLast) {
            for (size_t i {1}; i <= x.n_rows; ++i) {
                if (! seen.insert(x(x.n_rows - i)).second) {
                    // if duplicated, add index to vector res
                    res.push_back(i);
                }
            }
        } else {
            for (size_t i {0}; i < x.n_rows; ++i) {
                if (! seen.insert(x(i)).second) {
                    // if duplicated, add index to vector res
                    res.push_back(i);
                }
            }
        }
        return arma::conv_to<arma::uvec>::from(res);
    }
    // function that returns the indices of the first unique indices
    inline arma::uvec find_first_unique(const arma::vec& x,
                                        bool fromLast = false)
    {
        std::unordered_set<double> seen;
        std::vector<unsigned int> res;
        if (fromLast) {
            for (size_t i {1}; i <= x.n_rows; ++i) {
                if (seen.insert(x(x.n_rows - i)).second) {
                    // if duplicated, add index to vector res
                    res.push_back(i);
                }
            }
        } else {
            for (size_t i {0}; i < x.n_rows; ++i) {
                if (seen.insert(x(i)).second) {
                    // if duplicated, add index to vector res
                    res.push_back(i);
                }
            }
        }
        return arma::conv_to<arma::uvec>::from(res);
    }

    // set intersection for vector a and vector b
    // armadillo vector has just one template type parameter
    template <typename T, template <typename> class ARMA_VEC_TYPE>
    inline ARMA_VEC_TYPE<T> vec_intersection(const ARMA_VEC_TYPE<T>& a,
                                             const ARMA_VEC_TYPE<T>& b)
    {
        std::vector<T> res;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                              std::back_inserter(res));
        std::reverse(res.begin(), res.end());
        return arma::sort(arma::conv_to<ARMA_VEC_TYPE<T>>::from(res));
    }

    // set union for vector a and vector b
    template <typename T, template <typename> class ARMA_VEC_TYPE>
    inline ARMA_VEC_TYPE<T> vec_union(const ARMA_VEC_TYPE<T>& a,
                                      const ARMA_VEC_TYPE<T>& b)
    {
        std::vector<T> res;
        std::set_union(a.begin(), a.end(),
                       b.begin(), b.end(),
                       std::inserter(res, res.begin()));
        return arma::sort(arma::conv_to<ARMA_VEC_TYPE<T>>::from(res));
    }

    // cumulative sum in possibly reverse order
    inline arma::vec cum_sum(const arma::vec& x,
                             const bool reversely = false)
    {
        // if cumsum reversely
        if (reversely) {
            const unsigned long int n_x {x.n_rows};
            arma::vec res {arma::zeros(n_x)};
            double tmp {0.0};
            for (size_t i {1}; i <= n_x; ++i) {
                tmp += x[n_x - i];
                res[n_x - i] = tmp;
            }
            return res;
        }
        // otherwise, using arma::cumsum
        return arma::cumsum(x);
    }
    // cumulative sum of a matrix in possibly reverse order
    // dim = 0 for column-wise sum; dim = 1 for row-wise sum
    inline arma::mat cum_sum(const arma::mat& x,
                             const bool reversely = false,
                             unsigned int dim = 0)
    {
        // if cumsum reversely
        if (reversely) {
            if (dim == 0) {
                const unsigned long int n_x = x.n_rows;
                arma::mat tmp {x.row(n_x - 1)};
                arma::mat res {x};
                for (size_t i {2}; i <= n_x; ++i) {
                    tmp += x.row(n_x - i);
                    res.row(n_x - i) = tmp;
                }
                return res;
            } else if (dim == 1) {
                const unsigned long int n_x = x.n_cols;
                arma::mat tmp {x.col(n_x - 1)};
                arma::mat res {x};
                for (size_t i {2}; i <= n_x; ++i) {
                    tmp += x.col(n_x - i);
                    res.col(n_x - i) = tmp;
                }
                return res;
            } else {
                throw std::range_error("The 'dim' has to be either 0 or 1.");
            }
        }
        // otherwise, using arma::cumsum
        return arma::cumsum(x, dim);
    }
    // cumulative sum of a cube in possibly reverse order
    // along its third dimension
    inline arma::cube cum_sum(const arma::cube& x,
                              const bool reversely = false)
    {
        const unsigned long int n_x = x.n_slices;
        arma::cube res {x};
        // if cumsum reversely
        if (reversely) {
            arma::mat tmp {x.slice(n_x - 1)};
            for (size_t i {2}; i <= n_x; ++i) {
                tmp += x.slice(n_x - i);
                res.slice(n_x - i) = tmp;
            }
            return res;
        } else {
            arma::mat tmp {x.slice(0)};
            for (size_t i {1}; i < n_x; ++i) {
                tmp += x.slice(i);
                res.slice(i) = tmp;
            }
            return res;
        }
    }

    // aggregate sum of a vector based on same indices
    inline arma::vec aggregate_sum(const arma::vec& x,
                                   const arma::vec& indices,
                                   const bool simplify = true,
                                   const bool cumulative = false,
                                   const bool reversely = false)
    {
        const unsigned int n_x { x.n_elem };
        if (n_x != indices.n_elem) {
            throw std::logic_error(
                "The x and indices must have the same length."
                );
        }
        arma::vec uniInd { arma::unique(indices) };
        const unsigned int n_uniInd { uniInd.n_elem };

        // the x's having a same index are summed
        arma::vec sumVec { arma::zeros(n_uniInd) };

        // early exit if no need to aggregate for all unique indices
        bool is_all_unique { n_uniInd == n_x };
        if (is_all_unique) {
            arma::uvec sort_ind { arma::sort_index(indices) };
            sumVec = x.elem(sort_ind);
        } else {
            for (size_t i {0}; i < n_uniInd; ++i) {
                for (size_t j {0}; j < n_x; ++j) {
                    if (isAlmostEqual(uniInd[i], indices[j])) {
                        sumVec[i] += x[j];
                    }
                }
            }
        }
        if (cumulative) {
            sumVec = cum_sum(sumVec, reversely);
        }
        // if simplify the sum results to unique and sorted indices
        if (simplify || is_all_unique) {
            return sumVec;
        }
        // else
        arma::vec out {arma::zeros(n_x)};
        for (size_t i {0}; i < n_x; ++i) {
            for (size_t j {0}; j < n_uniInd; ++j) {
                if (isAlmostEqual(indices[i], uniInd[j])) {
                    out[i] = sumVec[j];
                    break;
                }
            }
        }
        return out;
    }
    // column-wise aggregrate sum for a matrix
    inline arma::mat aggregate_sum(const arma::mat& x,
                                   const arma::vec& indices,
                                   const bool simplify = true,
                                   const bool cumulative = false,
                                   const bool reversely = false)
    {
        const unsigned long int x_nrows {x.n_rows};
        if (x_nrows != indices.n_elem) {
            throw std::logic_error(
                "The number of rows of x must equal the length of indices."
                );
        }
        arma::vec uniInd {arma::unique(indices)};
        const unsigned long int n_uniInd {uniInd.n_rows};
        // the x's having a same index are summed
        arma::mat sumMat {arma::zeros(n_uniInd, x.n_cols)};
        for (size_t i {0}; i < n_uniInd; ++i) {
            for (size_t j {0}; j < x_nrows; ++j) {
                if (isAlmostEqual(uniInd[i], indices[j])) {
                    sumMat.row(i) += x.row(j);
                }
            }
        }
        if (cumulative) {
            sumMat = cum_sum(sumMat, reversely);
        }
        // if simplify the sum results to unique and sorted indices
        if (simplify) {
            return sumMat;
        }
        // else
        arma::mat out {arma::zeros(arma::size(x))};
        for (size_t i {0}; i < x_nrows; ++i) {
            for (size_t j {0}; j < n_uniInd; ++j) {
                if (isAlmostEqual(indices[i], uniInd[j])) {
                    out.row(i) = sumMat.row(j);
                    break;
                }
            }
        }
        return out;
    }
    // aggregrate sum for a cube over the third dimension
    inline arma::cube aggregate_sum(const arma::cube& x,
                                    const arma::vec& indices,
                                    const bool simplify = true,
                                    const bool cumulative = false,
                                    const bool reversely = false)
    {
        const unsigned long int x_nslices {x.n_slices};
        if (x_nslices != indices.n_elem) {
            throw std::logic_error(
                "The number of slices of x must equal the length of indices."
                );
        }
        arma::vec uniInd {arma::unique(indices)};
        const unsigned long int n_uniInd {uniInd.n_rows};
        // the x's having a same index are summed
        arma::cube sumCube {arma::zeros(x.n_rows, x.n_cols, n_uniInd)};
        for (size_t i {0}; i < n_uniInd; ++i) {
            for (size_t j {0}; j < x_nslices; ++j) {
                if (isAlmostEqual(uniInd[i], indices[j])) {
                    sumCube.slice(i) += x.slice(j);
                }
            }
        }
        if (cumulative) {
            sumCube = cum_sum(sumCube, reversely);
        }
        // if simplify the sum results to unique and sorted indices
        if (simplify) {
            return sumCube;
        }
        // else
        arma::cube out {arma::zeros(arma::size(x))};
        for (size_t i {0}; i < x_nslices; ++i) {
            for (size_t j {0}; j < n_uniInd; ++j) {
                if (isAlmostEqual(indices[i], uniInd[j])) {
                    out.slice(i) = sumCube.slice(j);
                    break;
                }
            }
        }
        return out;
    }


    // inline handy functions
    inline arma::vec mat2vec(const arma::mat& x) {
        return arma::conv_to<arma::vec>::from(x);
    }

    // function template for crossprod of two matrix-like objects
    template <typename T_matrix_like>
    inline arma::mat crossprod(const T_matrix_like& X,
                               const T_matrix_like& Y)
    {
        return X.t() * Y;
    }
    template <typename T_matrix_like>
    inline arma::mat crossprod(const T_matrix_like& X)
    {
        return X.t() * X;
    }
    // function template for tcrossprod of two matrix-like objects
    template <typename T_matrix_like>
    inline arma::mat tcrossprod(const T_matrix_like& X,
                                const T_matrix_like& Y)
    {
        return X * Y.t();
    }
    template <typename T_matrix_like>
    inline arma::mat tcrossprod(const T_matrix_like& X)
    {
        return X * X.t();
    }

    // function that computes L2-norm
    inline double l2_norm(const arma::vec& x)
    {
        return std::sqrt(arma::as_scalar(crossprod(x)));
    }
    inline double l2_norm(const arma::vec& x, const arma::vec& y)
    {
        return std::sqrt(arma::as_scalar(crossprod(x, y)));
    }

    // function computing relateive tolerance based on l2-norm
    inline double rel_l2_norm(const arma::vec& x_old, const arma::vec& x_new)
    {
        double denom { l2_norm(x_new + x_old) };
        if (isAlmostEqual(denom, 0)) {
            return 0;
        } else {
            return l2_norm(x_new - x_old) / denom;
        }
    }

    // function that computes L1-norm
    inline double l1_norm(const arma::vec& x)
    {
        return arma::sum(arma::abs(x));
    }

    // function computing relateive tolerance based on l1_norm
    inline double rel_l1_norm(const arma::vec& x_old, const arma::vec& x_new)
    {
        return l1_norm(x_new - x_old) / l1_norm(x_new + x_old);
    }


    // sign function
    inline double sign(const double& x)
    {
        if (x < 0) {
            return - 1.0;
        } else if (x > 0) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    // positive part
    template <typename T_scalar>
    inline T_scalar positive(T_scalar x)
    {
        if (x < 0) {
            return 0;
        } else {
            return x;
        }
    }

    // soft-thresholding operator
    inline double soft_threshold(const double& beta, const double& lambda)
    {
        return positive(std::abs(beta) - lambda) * sign(beta);
    }

    // convert uvec from logical comparison to uvec indices
    inline arma::uvec arma_which(const arma::uvec& x)
    {
        arma::uvec res { arma::zeros<arma::uvec>(arma::sum(x)) };
        for (size_t i {0}, j {0}; i < x.n_elem; ++i) {
            if (x[i] > 0) {
                res[j] = i;
                j++;
            }
        }
        return res;
    }

    // cumulative max in possibly reverse order
    inline arma::vec cum_max(const arma::vec& x,
                             const bool reversely = false)
    {
        arma::vec res { x };
        if (reversely) {
            for (size_t i {2}; i <= x.n_elem; ++i) {
                res(x.n_elem - i) = std::max(res(x.n_elem - i + 1),
                                             x(x.n_elem - i));
            }
        } else {
            for (size_t i {1}; i < x.n_elem; ++i) {
                res(i) = std::max(res(i - 1), x(i));
            }
        }
        return res;
    }
    inline arma::mat cum_max(const arma::mat& x,
                             const bool reversely = false)
    {
        arma::mat res { x };
        for (size_t i {0}; i < x.n_cols; ++i) {
            res.col(i) = cum_max(mat2vec(x.col(i)), reversely);
        }
        return res;
    }
    // cumulative min in possibly reverse order
    inline arma::vec cum_min(const arma::vec& x,
                             const bool reversely = false)
    {
        arma::vec res { x };
        if (reversely) {
            for (size_t i {2}; i <= x.n_elem; ++i) {
                res(x.n_elem - i) = std::min(res(x.n_elem - i + 1),
                                             x(x.n_elem - i));
            }
        } else {
            for (size_t i {1}; i < x.n_elem; ++i) {
                res(i) = std::min(res(i - 1), x(i));
            }
        }
        return res;
    }
    inline arma::mat cum_min(const arma::mat& x,
                             const bool reversely = false)
    {
        arma::mat res { x };
        for (size_t i {0}; i < x.n_cols; ++i) {
            res.col(i) = cum_min(mat2vec(x.col(i)), reversely);
        }
        return res;
    }

    // log of sum of exponentials
    inline double log_sum_exp(const arma::vec& x)
    {
        if (x.n_elem == 1) {
            return x(0);
        }
        double max_x { x.max() };
        double res { 0 };
        for (size_t i {0}; i < x.n_elem; ++i) {
            res += std::exp(x(i) - max_x);
        }
        return std::log(res) + max_x;
    }

    // repeat a number into a vector
    inline arma::vec rep_double(const double& x, const unsigned int& n) {
        return x * arma::ones(n);
    }
    inline arma::vec rep_each(const arma::vec& x, const unsigned int& n) {
        arma::vec res { arma::ones(n * x.n_elem) };
        for (size_t i {0}; i < res.n_elem; ++i) {
            size_t ind { i / n };
            res(i) = x(ind);
        }
        return res;
    }

    // step function
    inline arma::vec step_fun(const arma::vec& x,
                              const arma::vec& knots,
                              const arma::vec& height)
    {
        // create a map for fast comparison
        std::map<double, double> step_map;
        for (size_t i {0}; i < knots.n_elem; ++i) {
            step_map.insert(std::make_pair(knots(i), height(i + 1)));
        }
        arma::vec res { arma::zeros(x.n_elem) };
        std::map<double, double>::iterator it;
        for (size_t i {0}; i < x.n_elem; ++i) {
            it = step_map.upper_bound(x(i));
            if (it != step_map.begin()) {
                --it;
                res(i) = it->second;
            } else {
                res(i) = height(0);
            }
        }
        return res;
    }

    // quantile function
    // type 5 in quantile
    // reference: Hyndman and Fan (1996)
    inline double arma_quantile(const arma::vec& x, const double prob) {
        const double alpha { 0.5 };
        const unsigned int n { x.n_elem };
        if (prob < (1 - alpha) / n) {
            return x.min();
        }
        if (prob > (n - alpha) / n) {
            return x.max();
        }
        arma::vec inc_x { arma::sort(x) };
        int k { static_cast<int>(std::floor(n * prob + alpha)) };
        double pk { (k - alpha) / n };
        double w { (prob - pk) * n };
        return (1 - w) * inc_x(k - 1) + w * inc_x(k);
    }
    inline arma::vec arma_quantile(const arma::vec& x, const arma::vec& probs) {
        const double alpha { 0.5 };
        const unsigned int n { x.n_elem };
        arma::vec res { arma::zeros(probs.n_elem) };
        arma::vec inc_x { arma::sort(x) };
        for (size_t i {0}; i < probs.n_elem; ++i) {
            if (probs(i) < (1 - alpha) / n) {
                res(i) = x.min();
            } else if (probs(i) > (n - alpha) / n) {
                res(i) = x.max();
            } else {
                int k { static_cast<int>(std::floor(n * probs(i) + alpha)) };
                double pk { (k - alpha) / n };
                double w { (probs(i) - pk) * n };
                res(i) = (1 - w) * inc_x(k - 1) + w * inc_x(k);
            }
        }
        return res;
    }

    // convert arma vec type to Rcpp NumericVector
    template <typename T>
    Rcpp::NumericVector arma2rcppVec(const T& x) {
        return Rcpp::NumericVector(x.begin(), x.end());
    }

}

#endif
