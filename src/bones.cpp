#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]


bool AlmostRelativeEqual(double A, double B,
                         double MaxRelDiff = DBL_EPSILON)
{
  // Calculate the difference.
  double diff = fabs(A - B);
  A = fabs(A);
  B = fabs(B);
  // Find the largest
  double largest = (B > A) ? B : A;

  if (diff <= largest * MaxRelDiff)
    return true;
  return false;
}


// [[Rcpp::export]]
Rcpp::NumericVector aggregateSum(const Rcpp::NumericVector& x,
                                 const Rcpp::NumericVector& indices,
                                 bool simplify = true,
                                 bool cumulate = false,
                                 bool revcum = false)
{
  // naive function to sum x based on indices
  const unsigned long int n_x = x.size();

  Rcpp::NumericVector uniInd = Rcpp::sort_unique(indices);
  const unsigned long int n_uniInd = uniInd.size();

  // the x's having a same index are summed
  Rcpp::NumericVector sumVec(n_uniInd);
  for (size_t i = 0; i < n_uniInd; ++i) {
      for (size_t j = 0; j < n_x; ++j) {
        if (AlmostRelativeEqual(uniInd[i], indices[j])) {
          sumVec[i] += x[j];
        }
      }
  }
  // if `cumulate = true` for cumulated summation
  if (cumulate) {
    if (revcum) {
      long double tmp = sumVec[n_uniInd - 1];
      for (int i = n_uniInd - 2; i >= 0; --i) {
        tmp += sumVec[i];
        sumVec[i] = tmp;
      }
    } else {
      long double tmp = sumVec[0];
      for (size_t i = 1; i < n_uniInd; ++i) {
        tmp += sumVec[i];
        sumVec[i] = tmp;
      }
    }
  }
  // if simplify the sum results to unique and sorted indices
  if (simplify) {
    // add the names of indices? slows the function down a lot
    // if (addNames) {
    //   sumVec.attr("names") = Rcpp::as<Rcpp::CharacterVector>(uniInd);
    // }
    return sumVec;
  }

  // else
  Rcpp::NumericVector out(n_x);
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_uniInd; ++j) {
      if (AlmostRelativeEqual(indices[i], uniInd[j])) {
        out[i] = sumVec[j];
        break;
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector revcumsum(SEXP x)
{
  Rcpp::NumericVector xVec(x);
  const unsigned long int n_x = xVec.size();
  Rcpp::NumericVector res(n_x);
  long double tmp = 0.0;
  for (int i = n_x - 1; i >= 0; --i) {
    tmp += xVec[i];
    res[i] = tmp;
  }
  return res;
}
