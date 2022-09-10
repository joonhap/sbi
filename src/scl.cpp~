#include <Rcpp.h>
#include <iostream>
#include <chrono>

//' The SCL Distribution
//'
//' Quantile function, distribution function, and random generation for the SCL distribution. See Park and Won (2022) for information about the SCL distribution.
//'
//' @param p probability
//' @param n number of draws
//' @param M parameter to the SCL distribution
//' @param precision The requested level of precision for the outputs of qscl and pscl functions, in terms of the estimated standard deviation of the output. (Details: e.g., precision of 0.01 will output values with the standard deviation of approximately 0.01.)
//' @param lower logical; if TRUE, probabilities are P[X <= x], otherwise, P[X > x].
//' @param log_p logical; if TRUE, probabilities p are given as log(p).
//' @param force logical; if TRUE, the function will run regardless of how long it will take. If FALSE, the function will ask if you want to continue, stop, or give a new precision value whenever the expected run time is longer than 15 seconds. 
//' @return numeric vector of quantiles.
//' @examples
//' qscl(.99, 5, 2)
//' qscl(c(.01, .05, .95, .99), 10, 2.3)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qscl(Rcpp::NumericVector p, const double M, double precision, const bool lower = true, const bool log_p = false, const bool force = false) {
  int plen = p.size();
  if (log_p) {
    for (int i = 0; i < plen; i++) {
      p[i] = exp(p[i]);
    }
  }
  if (!lower) {
    for (int i = 0; i < plen; i++) {
      p[i] = 1.0 - p[i];
    }
  }
  
  const int nbloc = 50; // number of blocks
  int bsize = 200; // default block size
  std::vector<double> vec_scl(nbloc*bsize); // vector of SCL variates

  auto start = std::chrono::high_resolution_clock::now(); // start measuring time
  
  for (auto it = vec_scl.begin(); it != vec_scl.end(); it++) { // fill the random vector
    double z = R::rnorm(0.0, 1.0);
    double v = R::rchisq(M-1);
    *it = .5*(-z*z - v + M*log(v));
  }

  for (int i = 0; i < nbloc; i++) { // sort each block
    auto it = vec_scl.begin();
    std::sort(it, it + bsize);
    it += bsize;
  }

  auto stop = std::chrono::high_resolution_clock::now(); // stop the clock
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // duration in microseconds
  
  double max_sumsq = 0.0; // sum of squared deviations of quantiles computed from the blocks, and take maximum among plen such sumsq's.

  for (int i = 0; i < plen; i++) {
    double s = 0.0; // sum of quantiles
    double ss = 0.0; // sum of squared quantiles
    int r = round(p[i]*bsize) - 1;
    for (int j = 0; j < nbloc; j++) { // 100*p% quantiles for each block
      double q = vec_scl[j*bsize + r]; // quantile for the j-th block
      s += q;
      ss += q*q;
    }
    double sumsq = (ss - s*s/nbloc); // sum of squared deviations = (nbloc-1)*sample_variance
    if (sumsq > max_sumsq) { max_sumsq = sumsq; }
  }
  
  Rcpp::NumericVector quantiles(plen); // vector of quantiles
  double factor; // the factor by which the vector size increases (if necessary)
  int newbsize; // extended block size (if necessary)

  if (max_sumsq / (nbloc * (nbloc - 1))  < precision * precision) { // if the sample variance is sufficiently small, do ...
    std::sort(vec_scl.begin(), vec_scl.end()); // sort the entire vector
    for (int i = 0; i < plen; i++) {
      quantiles[i] = vec_scl[round(p[i]*bsize*nbloc) - 1];
    }
    return quantiles;
  } else { // else, increase the simulation length
    factor = max_sumsq / (nbloc * precision * precision * (nbloc - 1));
    double req_time = duration.count() / 1000.0 * factor; // estimated time for completion in seconds
    // std::cout << "The number of simulations will be increased by a factor of " << factor << '\n';
    if (req_time > 15.0 && (!force)) {
      do {
	std::cout << "This will take approximately " << round(duration.count() / 1000.0 * factor) << " seconds.\n";
	std::cout << "Do you want to continue? (If so, type 'y'.)\nIf not, you can enter a new precision (= the number of decimal points, but can be any positive real) or type 'n' to stop.\n";
	std::string response;
	std::cin >> response;
	if (response == "y" || response == "Y") {
	  break;
	} else if (response == "n" || response == "N") {
	  std::cout << "Stopping.\n";
	  return 0;
	} else {
	  precision = std::stod(response);
	  factor = max_sumsq / (nbloc * (nbloc - 1) * precision * precision);
	  continue;
	}
      } while (true);
    }
  }
  newbsize = round(factor * bsize);
  vec_scl.resize(nbloc*newbsize);
  
  for (int i = nbloc*bsize; i < nbloc*newbsize; i++) { // fill the extended block
    double z = R::rnorm(0.0, 1.0);
    double v = R::rchisq(M-1);
    vec_scl[i] = .5*(-z*z - v + M*log(v));
  }

  std::sort(vec_scl.begin(), vec_scl.end()); 
  for (int i = 0; i < plen; i++) {
    quantiles[i] = vec_scl[round(p[i]*newbsize*nbloc) - 1];
  }
  return quantiles;
}


// [[Rcpp::export]]
Rcpp::NumericVector pscl(Rcpp::NumericVector q, const double M, double precision, const bool lower = true, const bool log_p = false, const bool force = false) {
  int qlen = q.size();
  
  int vsize = 200; // starting size of the vector of random draws
  std::vector<double> vec_scl(vsize); // vector of SCL variates

  auto start = std::chrono::high_resolution_clock::now(); // start measuring time
  
  for (auto it = vec_scl.begin(); it != vec_scl.end(); it++) { // fill the random vector
    double z = R::rnorm(0.0, 1.0);
    double v = R::rchisq(M-1);
    *it = .5*(-z*z - v + M*log(v));
  }

  std::vector<int> counts(qlen); // number of random draws less than equal to q
  
  for (auto it = vec_scl.begin(); it != vec_scl.end(); it++) {
    for (int i = 0; i < qlen; i++) {
      counts[i] += (*it <= q[i]);
    }
  }
  
  auto stop = std::chrono::high_resolution_clock::now(); // stop the clock
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // duration in microseconds
  
  Rcpp::NumericVector probabilities(qlen); // vector of probabilities
  double factor; // the factor by which the vector size increases (if necessary)
  int newvsize; // extended vector size (if necessary)

  // revise from here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
  if (max_sumsq / (nbloc * (nbloc - 1))  < precision * precision) { // if the sample variance is sufficiently small, do ...
    std::sort(vec_scl.begin(), vec_scl.end()); // sort the entire vector
    for (int i = 0; i < plen; i++) {
      quantiles[i] = vec_scl[round(p[i]*bsize*nbloc) - 1];
    }
    return quantiles;
  } else { // else, increase the simulation length
    factor = max_sumsq / (nbloc * precision * precision * (nbloc - 1));
    double req_time = duration.count() / 1000.0 * factor; // estimated time for completion in seconds
    // std::cout << "The number of simulations will be increased by a factor of " << factor << '\n';
    if (req_time > 15.0 && (!force)) {
      do {
	std::cout << "This will take approximately " << round(duration.count() / 1000.0 * factor) << " seconds.\n";
	std::cout << "Do you want to continue? (If so, type 'y'.)\nIf not, you can enter a new precision (= the number of decimal points, but can be any positive real) or type 'n' to stop.\n";
	std::string response;
	std::cin >> response;
	if (response == "y" || response == "Y") {
	  break;
	} else if (response == "n" || response == "N") {
	  std::cout << "Stopping.\n";
	  return 0;
	} else {
	  precision = std::stod(response);
	  factor = max_sumsq / (nbloc * (nbloc - 1) * precision * precision);
	  continue;
	}
      } while (true);
    }
  }
  newbsize = round(factor * bsize);
  vec_scl.resize(nbloc*newbsize);
  
  for (int i = nbloc*bsize; i < nbloc*newbsize; i++) { // fill the extended block
    double z = R::rnorm(0.0, 1.0);
    double v = R::rchisq(M-1);
    vec_scl[i] = .5*(-z*z - v + M*log(v));
  }

  std::sort(vec_scl.begin(), vec_scl.end()); 
  for (int i = 0; i < plen; i++) {
    quantiles[i] = vec_scl[round(p[i]*newbsize*nbloc) - 1];
  }
  return quantiles;
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
