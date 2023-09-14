#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List PairwiseD(const NumericMatrix &cum, const NumericVector &b,
               const NumericVector &a){
  NumericMatrix xcum(cum);
  NumericVector xb(b);
  NumericVector xa(a);
  int n = xcum.ncol(), m = xa.size();
  NumericVector cum1, cum2, cum1p2, cum1m2;
  NumericVector cum1p22, cum1m22, xbb, delta, w;
  LogicalVector ind;
  NumericMatrix dmat(n, n * (m + 1));

  List ret;

  for (int i = 1; i < n; i++) {
    cum1 = xcum(_, i);
    for (int j = 0; j < i; j++) {

      cum2 = xcum(_, j);
      cum1p2 = cum1 + cum2;
      cum1m2 = Rcpp::abs(cum1 - cum2);
      ind = cum1p2 > 0;

      cum1p22 = cum1p2[ind];
      cum1m22 = cum1m2[ind];
      xbb = xb[ind];
      delta = cum1m22 / cum1p22;

      for (int k = 0; k < m; k++) {
        w = xbb * Rcpp::pow(cum1p22, xa[k]);
        dmat(i, j + k * n) = dmat(j, i + k * n) = Rcpp::sum(delta * w) / Rcpp::sum(w);
      }

      dmat(i, j + m * n) = dmat(j, i + m * n) = Rcpp::sum(Rcpp::ifelse(cum1p22 == cum1m22, xbb, 0)) / Rcpp::sum(xbb);

    }
  }
  ret["unifracs"] = dmat;

  return Rcpp::wrap(ret);
}

