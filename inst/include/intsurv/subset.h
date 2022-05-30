//
// intsurv: Integrative Survival Models
// Copyright (C) 2017-2022  Wenjie Wang <wang@wwenjie.org>
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

#ifndef INTSURV_SUBSET_H
#define INTSURV_SUBSET_H

#include <RcppArmadillo.h>
#include "CoxphCure.h"
#include "CoxphCureMar.h"

namespace intsurv {

    inline CoxphCure subset(const CoxphCure& object,
                            const arma::uvec& index)
    {
        CoxphCure out {
            object.surv_obj_.time_.elem(index),
            object.surv_obj_.event_.elem(index),
            object.surv_obj_.get_x(true, false).rows(index),
            object.cure_obj_.get_x(true, false).rows(index),
            object.control_,
            object.surv_obj_.control_,
            object.cure_obj_.control_
        };
        out.surv_obj_.set_offset(
            object.surv_obj_.control_.offset_.elem(index), false);
        const arma::uvec& out_ord { index(out.surv_obj_.ord_) };
        out.cure_obj_.set_offset(
            object.cure_obj_.control_.offset_.elem(out_ord)
            );
        return out;
    }

    inline CoxphCureMar subset(const CoxphCureMar& object,
                               const arma::uvec& index)
    {
        CoxphCureMar out {
            object.surv_obj_.time_.elem(index),
            object.get_event().elem(index),
            object.surv_obj_.get_x(true, false).rows(index),
            object.cure_obj_.get_x(true, false).rows(index),
            object.mar_obj_.get_x(true, false).rows(index),
            object.control_,
            object.surv_obj_.control_,
            object.cure_obj_.control_,
            object.mar_obj_.control_,
        };
        out.surv_obj_.set_offset(
            object.surv_obj_.control_.offset_.elem(index), false);
        const arma::uvec& out_ord { index(out.surv_obj_.ord_) };
        out.cure_obj_.set_offset(
            object.cure_obj_.control_.offset_.elem(out_ord)
            );
        out.mar_obj_.set_offset(
            object.mar_obj_.control_.offset_.elem(out_ord)
            );
        return out;
    }


}  // intsurv

#endif /* INTSURV_SUBSET_H */
