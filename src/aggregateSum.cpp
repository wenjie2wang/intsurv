#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector aggregateSum(Rcpp::NumericVector x,
                                 Rcpp::NumericVector indices,
                                 bool simplify = true,
                                 bool addNames = true)
{
  // naive function to sum x based on indices
  const long int n_x = x.size();

  Rcpp::NumericVector uniInd = Rcpp::sort_unique(indices);
  const long int n_uniInd = uniInd.size();

  // the x's having a same index are summed
  Rcpp::NumericVector sumVec(n_uniInd);
  for (size_t i = 0; i < n_uniInd; ++i) {
      for (size_t j = 0; j < n_x; ++j) {
        if (uniInd[i] == indices[j]) {
          sumVec[i] += x[j];
        }
      }
  }
  // if simplify the sum results to unique and sorted indices
  if (simplify) {
    // add the names of indices? slows the function down a lot
    if (addNames) {
      sumVec.attr("names") = Rcpp::as<Rcpp::CharacterVector>(uniInd);
    }
    return sumVec;
  }

  // else
  Rcpp::NumericVector out(n_x);
  for (size_t i = 0; i < n_x; ++i) {
    for (size_t j = 0; j < n_uniInd; ++j) {
      if (indices[i] == uniInd[j]) {
        out[i] = sumVec[j];
        break;
      }
    }
  }
  return out;
}
