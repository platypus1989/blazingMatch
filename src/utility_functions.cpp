// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' nearest score search c++ implementation.
//'
//' @param ctr_score matrix containing control group index and score
//' @param trt_score matrix containing treatment group index and score
//' @return matched treatment and control indices
//' @examples
//'
//' trt_score = data.frame(index=c(1:5), score=c(0.1, 0.03, 0.7, 0.2, 0.5))
//' ctr_score = data.frame(index=c(1:20), score=runif(20))
//' trt_score = trt_score[order(trt_score$score), ]
//' ctr_score = ctr_score[order(ctr_score$score), ]
//' nearestScore(as.matrix(ctr_score), as.matrix(trt_score))
//'
//' @export
// [[Rcpp::export]]
arma::mat nearestScore(arma::mat ctr_score,
                       arma::mat trt_score){
  int n0 = ctr_score.n_rows;
  int n1 = trt_score.n_rows;

  arma::mat results(n1, 2);

  int index0 = 0;
  int index1 = 0;

  while (index0 < n0 && index1 < n1) {
    while (index0 + 1 < n0 && ctr_score(index0+1, 1) < trt_score(index1, 1)) {
      index0++;
    }
    if (index0 + 1 == n0) {
      for (int i=index1; i< n1; i++) {
        results(i, 0) = ctr_score(index0, 0);
        results(i, 1) = trt_score(i, 0);
      }
      return results;
    }
    if (ctr_score(index0+1, 1) - trt_score(index1, 1) >= trt_score(index1, 1) - ctr_score(index0, 1)) {
      results(index1, 0) = ctr_score(index0, 0);
      results(index1, 1) = trt_score(index1, 0);
    } else {
      results(index1, 0) = ctr_score(index0+1, 0);
      results(index1, 1) = trt_score(index1, 0);
      index0++;
    }
    index1++;
  }
  return results;
}
