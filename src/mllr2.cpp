#include <Rcpp.h>
#include <iostream>
#include <chrono>

//' The MLLR_2 distribution
//' 
//' Quantile function, distribution function, and random generation for the MLLR_2 distribution family. See Park (2023) for information about the MLLR distributions.
//'
//' @name MLLR2
//' @param p vector of probabilities
//' @param q vector of quantiles
//' @param n number of draws
//' @param M the first parameter for the MLLR_2 distribution
//' @param k the second parameter for the MLLR_2 distribution
//' @param num_error_size The requested size of numerical error for the outputs of qmllr2 and pmllr2 functions, in terms of the estimated standard deviation of the output. For example num_error_size of 0.01 will output values with the standard deviation of approximately equal to 0.01.
//' @param lower logical; if TRUE, probabilities are P[X <= x], otherwise, P[X > x].
//' @param log_p logical; if TRUE, probabilities p are given as log(p).
//' @param force logical; if TRUE, the function will run regardless of how long it will take. If FALSE, the function will ask if you want to continue, stop, or give a new num_error_size value whenever the expected run time is longer than 15 seconds. 
//' @return a list consisting of the numeric vector of quantiles and the num_error_size (numeric) used.
//' @examples
//' qmllr2(.99, 5, 2, 1.2, 0.23)
//' qmllr2(c(.01, .05, .95, .99), 10, 2.3, 3, 10)
//' qmllr2(c(.01, .05, .95, .99), 10, 2.3, 3, 10, num_error_size=0.01, lower=TRUE)
//' pmllr2(c(-8.3, -5.9), 8, 1, 3.8, 4.1)
//' pmllr2(c(-8.3, -5.9), 8 ,1, 3.8, 4.1, force=TRUE)
//' rmllr2(10, 7, 2, 91, 1.1)
//' @export
// [[Rcpp::export]]
Rcpp::List qmllr2(Rcpp::NumericVector p, const double M, const double k, const double nu, const double s0, double num_error_size = 0.01, const bool lower = true, const bool log_p = false, const bool force = false) {
  if (k <= 0) { 
    std::cout << "k should be positive." << std::endl;
    exit(1);
  }
  if (M <= k) { 
    std::cout << "M should be greater than k." << std::endl;
    exit(1);
  }
  if (nu <= 0) { 
    std::cout << "nu should be positive." << std::endl;
    exit(1);
  }
  if (s0 <= 0) { 
    std::cout << "s0 should be positive." << std::endl;
    exit(1);
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

  // draw from the MLLR_2 distribution in blocks, compute the quantile for each block, and check the standard deviation of the quantiles.
  const int nbloc = 50; // number of blocks
  int bsize = 200; // default block size
  std::vector<double> vec_mllr2(nbloc*bsize); // vector of MLLR_2 variates

  auto start = std::chrono::high_resolution_clock::now(); // start measuring time
  
  for (auto it = vec_mllr2.begin(); it != vec_mllr2.end(); it++) { // fill the random vector
    double z1 = R::rnorm(0,1);
    double x2 = R::rchisq(M-k);
    double z1s = z1 - s0/2/sqrt(nu);
    double sstsq = 2.0*nu*M*(sqrt(s0*s0/(nu*M*M)*(z1s*z1s+x2)+1)-1);
    double z1ss = z1s + sstsq/2/sqrt(nu)/s0;
    *it = -s0*s0/2/sstsq*(z1ss*z1ss+x2) + M/2.0*log(s0*s0*x2/(M*sstsq)) + M/2.0;
  }

  for (int i = 0; i < nbloc; i++) { // sort each block
    auto it = vec_mllr2.begin();
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
      double q = vec_mllr2[j*bsize + r]; // quantile for the j-th block
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
    std::sort(vec_mllr2.begin(), vec_mllr2.end()); // sort the entire vector
    for (int i = 0; i < plen; i++) {
      quantiles[i] = vec_mllr2[round(p[i]*bsize*nbloc) - 1];
    }
    return Rcpp::List::create(Rcpp::Named("quantiles") = quantiles, Rcpp::Named("numerical_error_size") = num_error_size);
  } else { // else, increase the simulation length
    factor = max_sumsq / (nbloc * num_error_size * num_error_size * (nbloc - 1));
    double req_time = duration.count() / 1000.0 * factor; // estimated time for completion in seconds

    if (req_time > 15.0 && (!force)) {
      do {
	std::cout << "Computing quantile values for the MLLR_2 distribution (" << M << "," << k << "," << nu << "," << s0  << ") with approximate size of numerical error " << num_error_size << ".\n";
	std::cout << "This will take approximately " << round(duration.count() / 1000.0 * factor) << " seconds.\n";
	std::cout << "Do you want to continue? (If so, type 'y'.)\nIf not, you can enter a new target approximate numerical error size (e.g., 0.03) or type 'n' to stop.\n";
	std::string response;
	std::cin >> response;
	if (response == "y" || response == "Y") {
	  break;
	} else if (response == "n" || response == "N") {
	  std::cout << "Stopping.\n";
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
  vec_mllr2.resize(nbloc*newbsize);
  
  for (int i = nbloc*bsize; i < nbloc*newbsize; i++) { // fill the extended block
    double z1 = R::rnorm(0,1);
    double x2 = R::rchisq(M-k);
    double z1s = z1 - s0/2/sqrt(nu);
    double sstsq = 2.0*nu*M*(sqrt(s0*s0/(nu*M*M)*(z1s*z1s+x2)+1)-1);
    double z1ss = z1s + sstsq/2/sqrt(nu)/s0;
    vec_mllr2[i] = -s0*s0/2/sstsq*(z1ss*z1ss+x2) + M/2.0*log(s0*s0*x2/(M*sstsq)) + M/2.0;
  }

  std::sort(vec_mllr2.begin(), vec_mllr2.end()); 
  for (int i = 0; i < plen; i++) {
    quantiles[i] = vec_mllr2[round(p[i]*newbsize*nbloc) - 1];
  }

  return Rcpp::List::create(Rcpp::Named("quantiles") = quantiles, Rcpp::Named("numerical_error_size") = num_error_size);
}




//' @rdname MLLR2
// [[Rcpp::export]]
Rcpp::List pmllr2(Rcpp::NumericVector q, const double M, const double k, const double nu, const double s0, double num_error_size = 0.01, const bool lower = true, const bool log_p = false, const bool force = false) {
  if (k <= 0) { 
    std::cout << "k should be positive." << std::endl;
    exit(1);
  }
  if (M <= k) { 
    std::cout << "M should be greater than k." << std::endl;
    exit(1);
  }
  if (nu <= 0) { 
    std::cout << "nu should be positive." << std::endl;
    exit(1);
  }
  if (s0 <= 0) { 
    std::cout << "s0 should be positive." << std::endl;
    exit(1);
  }

  int qlen = q.size();
  
  int vsize = std::max(2000, int(1/num_error_size)); // starting size of the vector of random draws
  std::vector<int> counts(qlen); // number of random draws less than equal to q
  
  auto start = std::chrono::high_resolution_clock::now(); // start measuring time
  
  for (int i = 0; i < vsize; i++) { // generate MLLR_2 random variates
    double z1 = R::rnorm(0,1);
    double x2 = R::rchisq(M-k);
    double z1s = z1 - s0/2/sqrt(nu);
    double sstsq = 2.0*nu*M*(sqrt(s0*s0/(nu*M*M)*(z1s*z1s+x2)+1)-1);
    double z1ss = z1s + sstsq/2/sqrt(nu)/s0;
    double draw = -s0*s0/2/sstsq*(z1ss*z1ss+x2) + M/2.0*log(s0*s0*x2/(M*sstsq)) + M/2.0;
    for (int j = 0; j < qlen; j++) {
      counts[j] += (draw <= q[j]);
    }
  }

  auto stop = std::chrono::high_resolution_clock::now(); // stop the clock
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // duration in microseconds
  
  Rcpp::NumericVector probs(qlen); // vector of probabilities

  double max_var = 0.0;
  for (int i = 0; i < qlen; i++) {
    if (double(counts[i]) * (vsize - counts[i]) > max_var) { max_var = double(counts[i]) * (vsize - counts[i]); }
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
	std::cout << "Computing the cdf for the MLLR_2 distribution (" << M << "," << k << ") with approximate size of numerical error " << num_error_size << ".\n";
	std::cout << "This will take approximately " << round(duration.count() / 1000.0 * factor) << " seconds.\n";
	std::cout << "Do you want to continue? (If so, type 'y'.)\nIf not, you can enter a new approximate size for numerical error (e.g., 0.03) or type 'n' to stop.\n";
	std::string response;
	std::cin >> response;
	if (response == "y" || response == "Y") {
	  break;
	} else if (response == "n" || response == "N") {
	  std::cout << "Stopping.\n";
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
    double z1 = R::rnorm(0,1);
    double x2 = R::rchisq(M-k);
    double z1s = z1 - s0/2/sqrt(nu);
    double sstsq = 2.0*nu*M*(sqrt(s0*s0/(nu*M*M)*(z1s*z1s+x2)+1)-1);
    double z1ss = z1s + sstsq/2/sqrt(nu)/s0;
    double draw = -s0*s0/2.0/sstsq*(z1ss*z1ss+x2) + M/2.0*log(s0*s0*x2/M/sstsq) + M/2.0;

    for (int j = 0; j < qlen; j++) {
      counts[j] += (draw <= q[j]);
    }
  }

  for (int j = 0; j < qlen; j++) {
    probs[j] = counts[j] / (newvsize + 0.0);
  }

  return Rcpp::List::create(Rcpp::Named("probs") = probs, Rcpp::Named("numerical_error_size") = num_error_size);
}




//' @rdname MLLR2
// [[Rcpp::export]]
Rcpp::NumericVector rmllr2(const int n, const double M, const double k, const double nu, const double s0) {
  if (k <= 0) { 
    std::cout << "k should be positive." << std::endl;
    exit(1);
  }
  if (M <= k) { 
    std::cout << "M should be greater than k." << std::endl;
    exit(1);
  }
  if (nu <= 0) { 
    std::cout << "nu should be positive." << std::endl;
    exit(1);
  }
  if (s0 <= 0) { 
    std::cout << "s0 should be positive." << std::endl;
    exit(1);
  }
  
  Rcpp::NumericVector out(n); // output

  for (int i = 0; i < n; i++) {
    double z1 = R::rnorm(0,1);
    double x2 = R::rchisq(M-k);
    double z1s = z1 - s0/2/sqrt(nu);
    double sstsq = 2.0*nu*M*(sqrt(s0*s0/(nu*M*M)*(z1s*z1s+x2)+1)-1);
    double z1ss = z1s + sstsq/2/sqrt(nu)/s0;
    out[i] = -s0*s0/2.0/sstsq*(z1ss*z1ss+x2) + M/2.0*log(s0*s0*x2/M/sstsq) + M/2.0;
  }

  return out;
}
