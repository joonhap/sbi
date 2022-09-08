#include <Rcpp.h>

//' The SCL Distribution
//'
//' Quantile function, distribution function, and random generation for the SCL distribution. See Park and Won (2022) for information about the SCL distribution.
//'
//' @param p probability
//' @param M parameter to the SCL distribution
//' @param precision The requested level of precision for the outputs of qscl and pscl functions, as the number of decimal places. See Details for more information. (Details: precision can be any positive real, and checks if the standard error of the output is not greater than 10^{-precision}.)
//' @return b b
//' @examples
//' asdf5
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qscl(Rcpp::NumericVector p, double M, double precision, bool lower = true, bool log = false) {
  const int nbloc = 50; // number of blocks
  int bsize = 500; // default block size
  std::vector<double> vec_scl(nbloc*bsize); // vector of SCL variates

  for (auto it : vec_scl) { // fill the random vector
    double z = R::rnorm(0.0, 1.0);
    double v = R::rchisq(M-1);
    *it = .5*(-z*z - v + M*log(v));
  }

  for (int i = 0; i < nbloc; i++) { // sort each block
    auto it = vec_scl.begin();
    std::sort(it, it + bsize);
    it += bsize;
  }

  std::vector<double> quantiles(nbloc); // vector of quantiles
  
}

// [[Rcpp::export]]
Rcpp::NumericVector rscl(int n, double M) {
  Rcpp::NumericVector out(n); // output

  for (int i = 0; i < n; i++) {
    double z = R::rnorm(0.0, 1.0);
    double v = R::rchisq(M-1);
    out[i] = .5*(-z*z - v + M*log(v));
  }

  return out;
}
