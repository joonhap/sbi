#include <Rcpp.h>
#include <iostream>
#include <chrono>

//' The SCL distribution
//'
//' Quantile function, distribution function, and random generation for the SCL distribution family. See Park (2025) for information about the SCL distributions.
//'
//' @name SCL
//' @param p vector of probabilities
//' @param q vector of quantiles
//' @param n number of draws
//' @param M the first parameter for the SCL distributions
//' @param k the second parameter for the SCL distribution
//' @param num_error_size The requested size of numerical error for the outputs of qscl and pscl functions, in terms of the estimated standard deviation of the output. For example num_error_size of 0.01 will output values with the standard deviation of approximately equal to 0.01.
//' @param lower logical; if TRUE, probabilities are P(X <= x), otherwise, P(X > x).
//' @param log_p logical; if TRUE, probabilities p are given as log(p).
//' @param force logical; if TRUE, the function will run regardless of how long it will take. If FALSE, the function will ask if you want to continue, stop, or give a new num_error_size value whenever the expected run time is longer than 15 seconds. 
//' @return a list consisting of the numeric vector of quantiles and the num_error_size (numeric) used.
//' @export
// [[Rcpp::export]]
Rcpp::List qscl(Rcpp::NumericVector p, const double M, const double k, double num_error_size = 0.01, const bool lower = true, const bool log_p = false, const bool force = false) {
  if (k <= 0) { 
    Rcpp::stop("k should be positive.\n");
  }
  if (M <= k) { 
    Rcpp::stop("M should be greater than k.\n");
  }
  
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

  // draw from the SCL distribution in blocks, compute the quantile for each block, and check the standard deviation of the quantiles.
  const int nbloc = 50; // number of blocks
  int bsize = 200; // default block size
  std::vector<double> vec_scl(nbloc*bsize); // vector of SCL variates

  auto start = std::chrono::high_resolution_clock::now(); // start measuring time
  
  for (auto it = vec_scl.begin(); it != vec_scl.end(); it++) { // fill the random vector
    double x1 = R::rchisq(k);
    double x2 = R::rchisq(M-k);
    *it = .5*(-x1 - x2 + M*log(x2/M) + M);
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

  if (max_sumsq / (nbloc * (nbloc - 1))  < num_error_size * num_error_size) { // if the sample variance is sufficiently small, do ...
    std::sort(vec_scl.begin(), vec_scl.end()); // sort the entire vector
    for (int i = 0; i < plen; i++) {
      quantiles[i] = vec_scl[round(p[i]*bsize*nbloc) - 1];
    }
    return Rcpp::List::create(Rcpp::Named("quantiles") = quantiles, Rcpp::Named("numerical_error_size") = num_error_size);
  } else { // else, increase the simulation length
    factor = max_sumsq / (nbloc * num_error_size * num_error_size * (nbloc - 1));
    double req_time = duration.count() / 1000.0 * factor; // estimated time for completion in seconds

    if (req_time > 15.0 && (!force)) {
      do {
	Rcpp::Rcout << "Computing quantile values for the SCL distribution (" << M << "," << k << ") with approximate size of numerical error " << num_error_size << ".\n";
	Rcpp::Rcout << "This will take approximately " << round(duration.count() / 1000.0 * factor) << " seconds.\n";
	Rcpp::Rcout << "Do you want to continue? (If so, type 'y'.)\nIf not, you can enter a new approximate numerical error size (e.g., 0.03) or type 'n' to stop.\n";
	std::string response;
	std::cin >> response;
	if (response == "y" || response == "Y") {
	  break;
	} else if (response == "n" || response == "N") {
	  Rcpp::Rcout << "Stopping.\n";
	  return 0;
	} else {
	  num_error_size = std::stod(response);
	  factor = max_sumsq / (nbloc * (nbloc - 1) * num_error_size * num_error_size);
	  continue;
	}
      } while (duration.count() / 1000.0 * factor > 15.0);
    }
  }
  newbsize = round(factor * bsize);
  vec_scl.resize(nbloc*newbsize);

  for (int i = nbloc*bsize; i < nbloc*newbsize; i++) { // fill the extended block
    double x1 = R::rchisq(k);
    double x2 = R::rchisq(M-k);
    vec_scl[i] = .5*(-x1 - x2 + M*log(x2/M) + M);
  }

  std::sort(vec_scl.begin(), vec_scl.end()); 
  for (int i = 0; i < plen; i++) {
    quantiles[i] = vec_scl[round(p[i]*newbsize*nbloc) - 1];
  }

  return Rcpp::List::create(Rcpp::Named("quantiles") = quantiles, Rcpp::Named("numerical_error_size") = num_error_size);
}




//' @rdname SCL
//' @export
// [[Rcpp::export]]
Rcpp::List pscl(Rcpp::NumericVector q, const double M, const double k, double num_error_size = 0.01, const bool lower = true, const bool log_p = false, const bool force = false) {
  if (k <= 0) { 
    Rcpp::stop("k should be positive.\n");
  }
  if (M <= k) { 
    Rcpp::stop("M should be greater than k.\n");
  }
  
  int qlen = q.size();
  
  int vsize = std::max(2000, int(1/num_error_size)); // starting size of the vector of random draws
  std::vector<int> counts(qlen); // number of random draws less than equal to q
  
  auto start = std::chrono::high_resolution_clock::now(); // start measuring time
  
  for (int i = 0; i < vsize; i++) { // generate SCL random variates
    double x1 = R::rchisq(k);
    double x2 = R::rchisq(M-k);
    double draw = .5*(-x1 - x2 + M*log(x2/M) + M);

    for (int j = 0; j < qlen; j++) {
      counts[j] += (draw <= q[j]);
    }
  }

  auto stop = std::chrono::high_resolution_clock::now(); // stop the clock
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // duration in microseconds
  
  Rcpp::NumericVector probs(qlen); // vector of probabilities

  double max_var = 0.0;
  for (int i = 0; i < qlen; i++) {
    if (double(counts[i]) * (vsize - counts[i]) > max_var) {
      max_var = double(counts[i]) * (vsize - counts[i]);
    }
  }

  max_var = std::max(max_var, double(vsize)) / (double(vsize)*vsize*vsize);
  double factor = max_var / (num_error_size * num_error_size); // the factor by which the vector size increases (if necessary) 

  if (factor < 1.0) { // if all sample std.dev is smaller than num_error_size
    for (int i = 0; i < qlen; i++) {
      probs[i] = counts[i] / (vsize + 0.0);
    }
    return Rcpp::List::create(Rcpp::Named("probs") = probs, Rcpp::Named("numerical_error_size") = num_error_size);
  } else { // else, increase the number of draws
    double req_time = duration.count() / 1000.0 * factor; // estimated time for completion in seconds

    if (req_time > 15.0 && (!force)) {
      do {
	Rcpp::Rcout << "Computing the cdf for the SCL distribution (" << M << "," << k << ") with approximate size of numerical error " << num_error_size << ".\n";
	Rcpp::Rcout << "This will take approximately " << round(duration.count() / 1000.0 * factor) << " seconds.\n";
	Rcpp::Rcout << "Do you want to continue? (If so, type 'y'.)\nIf not, you can enter a new approximate numerical error size (e.g., 0.03) or type 'n' to stop.\n";
	std::string response;
	std::cin >> response;
	if (response == "y" || response == "Y") {
	  break;
	} else if (response == "n" || response == "N") {
	  Rcpp::Rcout << "Stopping.\n";
	  return 0;
	} else {
	  double num_error_size = std::stod(response);
	  factor = max_var / (num_error_size * num_error_size);
	  continue;
	}
      } while (duration.count() / 1000.0 * factor > 15.0);
    }
  }
  int newvsize = round(factor * vsize); // extended vector size (if necessary)
  
  for (int i = vsize; i < newvsize; i++) { // fill the extended vector with random draws
    double x1 = R::rchisq(k);
    double x2 = R::rchisq(M-k);
    double draw = .5*(-x1 - x2 + M*log(x2/M) + M);

    for (int j = 0; j < qlen; j++) {
      counts[j] += (draw <= q[j]);
    }
  }

  for (int j = 0; j < qlen; j++) {
    probs[j] = counts[j] / (newvsize + 0.0);
  }

  return Rcpp::List::create(Rcpp::Named("probs") = probs, Rcpp::Named("numerical_error_size") = num_error_size);
}




//' @rdname SCL
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rscl(const int n, const double M, const double k) {
  if (k <= 0) { 
    Rcpp::stop("k should be positive.\n");
  }
  if (M <= k) { 
    Rcpp::stop("M should be greater than k.\n");
  }
  
  Rcpp::NumericVector out(n); // output

  for (int i = 0; i < n; i++) {
    double x1 = R::rchisq(k);
    double x2 = R::rchisq(M-k);
    out[i] = .5*(-x1 - x2 + M*log(x2/M) + M);
  }

  return out;
}
