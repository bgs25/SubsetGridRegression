// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;



NumericVector colMax(NumericMatrix mat, bool na_rm=false) {
  NumericVector out = mat(0, _);
  if (na_rm) {
    for (int j=0; j<out.size(); j++) {
      if (NumericVector::is_na(out[j])) out[j] = R_NegInf;
    }
    for (int j=0; j<mat.ncol(); j++) {
      for (int i=1; i<mat.nrow(); i++) {
        if (!NumericVector::is_na(mat(i, j))) {
          if (mat(i, j) > out[j]) out[j] = mat(i, j);
        }
      }
    }
  }
  else {
    for (int j=0; j<mat.ncol(); j++) {
      for (int i=1; i<mat.nrow(); i++) {
        if (mat(i, j) > out[j]) out[j] = mat(i, j);
      }
    }
  }
  return out;
}


// Functions to compute the Lasso, square root Lasso and and Lasso with multiple responses
// Uses a vector of starting values and initialisers for the active est, strong set and etc.
// Always uses covariance updates
// Templates used to allow functions to be applied to vectors or columns of matrices without
// copying. Templated versions typically end _c with the exported equivalents ending _R
// Typically templated parts must be either a matrix column or vector

inline double soft_thresh(double u, double t) {
  if (u > 0) {
    if (u > t) {
      return u-t;
    }
  } else {
    if (u < -t) {
      return u+t;
    }
  }
  return 0;
}

// Computes ||Y - X\beta||_2^2
template<typename T1, typename T2>
double comp_RSS(double const null_dev, T1 const Xty, NumericMatrix const XtX,
                T2 const beta,
                IntegerVector const A, int const A_size) {
  double RSS=0;
  double cross_prod=0;
  double Xbeta_sq=0;
  for (int k=0; k<A_size; k++) {
    cross_prod += beta[A[k]] * Xty[A[k]];
    RSS += XtX(A[k], A[k]) * beta[A[k]] * beta[A[k]];
    for (int j=0; j<k; j++) {
      Xbeta_sq += XtX(A[j], A[k]) * beta[A[k]] * beta[A[j]];
    }
  }
  RSS += 2*(Xbeta_sq - cross_prod) + null_dev;
  return RSS;
}

// Rescales beta by a factor constant across components according to the least squares regression of the response
// on to the "direction" given by beta. This helps to debias beta and improve prediction error
template<typename T1, typename T2>
double rescale_fit(double const null_dev, T1 const Xty, NumericMatrix const XtX,
                   T2 beta, IntegerVector const A, int const A_size) {
  double Xbeta_sq_pure=0;
  double cross_prod=0;
  double Xbeta_sq_cross=0;
  for (int k=0; k<A_size; k++) {
    cross_prod += beta[A[k]] * Xty[A[k]];
    Xbeta_sq_pure += XtX(A[k], A[k]) * beta[A[k]] * beta[A[k]]; // pure quadratic terms
    for (int j=0; j<k; j++) {
      Xbeta_sq_cross += XtX(A[j], A[k]) * beta[A[k]] * beta[A[j]];
    }
  }
  double Xbeta_sq = Xbeta_sq_pure + 2*Xbeta_sq_cross;
  double scale_fac = cross_prod / Xbeta_sq;
  // rescale beta
  for (int k=0; k<A_size; k++) {
    beta[A[k]] *= scale_fac;
  }
  double RSS = Xbeta_sq - 2*cross_prod + null_dev;
  return RSS;
}


double rescale_fit_R(double const null_dev, NumericVector Xty, NumericMatrix const XtX,
                     NumericVector beta, IntegerVector const A, int const A_size) {
  return rescale_fit(null_dev, Xty, XtX, beta, A, A_size);
}


int rescale_fit_loop(NumericVector const null_dev, NumericMatrix XtY,
                     NumericMatrix const XtX, NumericMatrix beta) {
  int p=beta.nrow();
  for (int k=0; k<null_dev.size(); k++) {
    // Find active set
    IntegerVector A(p);
    int A_size=0;
    for (int j=0; j<p; j++) {
      if (beta(j, k) != 0) {
        A[A_size] = j;
        A_size++;
      }
    }
    NumericMatrix::Column betacol = beta( _, k);
    NumericMatrix::Column Xty = XtY(_, k);
    
    rescale_fit(null_dev[k], Xty, XtX, betacol, A, A_size);
  }
  return 0;
}


// not used yet
template<typename T1, typename T2, typename T3>
int comp_resid(T1 resid_out, T2 const y, NumericMatrix const x, T3 const beta,
               NumericVector const A, NumericVector const A_size) {
  for (int k=0; k<A_size[0]; k++) {
    int A_k = A[k];
    for (int i=0; i<resid_out.size(); i++) {
      resid_out[i] += beta[A_k] * x(i, A_k);
    }
  }
  return 0;
}

// thresh should be 1e-7 * null_dev 
// X and y centred
// Uses covariance updates
// Updates beta
// active set A is zero-indexed
template<typename T1, typename T2>
int beta_active(T1 beta, T2 const Xty, NumericMatrix const XtX,
                IntegerVector const A, int const A_size, double const thresh, int const maxit,
                double lam) {
  // l will be 0 indexed
  // Could maintain a gradient vector
  // Rcout << "start beta_active" << std::endl;
  double rel_err;
  int iter=0;
  do {
    rel_err=0;
    
    for (int k=0; k<A_size; k++) {
      int A_k=A[k];
      double X_ktR = Xty[A_k];
      for (int j=0; j<A_size; j++) {
        X_ktR -= XtX(A[j], A_k) * beta[A[j]];
      }
      
      // Update beta
      double new_coef = soft_thresh(X_ktR + XtX(A_k, A_k) * beta[A_k], lam) / XtX(A_k, A_k);
      double diff = new_coef - beta[A_k];
      double abs_diff = std::abs(beta[A_k]) - std::abs(new_coef);
      beta[A_k] = new_coef;
      
      
      
      // Update rel_err
      double obj_change = diff*(X_ktR - diff*XtX(A_k, A_k)/2) + lam*abs_diff;
      
      rel_err = std::max(obj_change, rel_err);
    }
    iter++;
  } while (rel_err > thresh && iter < maxit);
  
  // Rcout << "end beta_active" << std::endl;
  return 0;
}

// 0 means no changes
// 1 or more means changes were made so beta_active needs to be run once more
// S, A  and R will always be disjoint (R is the "rest" of the predictors)
// Implements KKT check at the bottom of pg 20 of
// https://statweb.stanford.edu/~tibs/ftp/strong.pdf

template<typename T1, typename T2>
int KKT_check(IntegerVector A, IntegerVector A_size, IntegerVector S, IntegerVector S_size,
              IntegerVector R, IntegerVector R_size, T1 const beta,
              T2 const Xty, NumericMatrix const XtX,
              double const lam, double const lam_next) {
  int out=0;
  //REMOVE//////21SEP Rcout << "KKT check" << std::endl;
  // Check for violations in strong set S
  for (int k=0; k<S_size[0]; k++) {
    double X_ktR = Xty[S[k]];
    for (int j=0; j<A_size[0]; j++) {
      X_ktR -= XtX(A[j], S[k]) * beta[A[j]];
    }
    // Add violations to the strong set S, and remove these from S
    if (std::abs(X_ktR) > lam) {
      A[A_size[0]++] = S[k]; // Add elem to A
      S[k] = S[--S_size[0]]; // Remove kth elem of S
      out=1;
    }
  }
  
  if (out==1) {
    //REMOVE//////21SEP Rcout << "end KKT check; found violator in strong set" << std::endl;
    return out;
  }
  
  double strong_thresh = 2*lam_next - lam;
  // Check for violations among all predictors (there are no violations in S or A)
  for (int k=0; k<R_size[0]; k++) {
    double X_ktR = Xty[R[k]];
    for (int j=0; j<A_size[0]; j++) {
      X_ktR -= XtX(A[j], R[k]) * beta[A[j]];
    }
    if (std::abs(X_ktR) > lam) {
      A[A_size[0]++] = R[k]; // Add elem to A
      R[k] = R[--R_size[0]]; // Remove kth elem of R
      out=2;
    } else if (std::abs(X_ktR) >= strong_thresh) {
      S[S_size[0]++] = R[k]; // Add elem to S
      R[k] = R[--R_size[0]]; // Remove kth elem of R
    }
  }
  
  return out;
}



template<typename T1, typename T2, typename T3> // look up templates BGS
int KKT_check_m (T3 A_add_vec, int& A_add_size_num,
                 IntegerVector A, IntegerVector A_size, IntegerVector S, IntegerVector S_size,
                 IntegerVector R, IntegerVector R_size, T1 const beta,
                 T2 const Xty, NumericMatrix const XtX,
                 double const lam, double const lam_next) {
  int out=0;
  // Check for violations in strong set S
  for (int k=0; k<S_size[0]; k++) {
    double X_ktR = Xty[S[k]];
    for (int j=0; j<A_size[0]; j++) {
      X_ktR -= XtX(A[j], S[k]) * beta[A[j]];
    }
    // Add violations to the strong set S, and remove these from S
    if (std::abs(X_ktR) > lam) {
      A[A_size[0]++] = S[k]; // Add elem to A
      A_add_vec[ A_add_size_num++ ]  = S[ k ]; // Need to pass l, A_add and A_add_size!
      S[k] = S[ --S_size[0] ]; // Remove kth elem of S
      out=1;
    }
  }
  
  if (out==1) {
    return out;
  }
  
  double strong_thresh = 2*lam_next - lam;
  // Check for violations among all predictors (there are no violations in S or A)
  
  for (int k=0; k<R_size[0]; k++) {
    double X_ktR = Xty[R[k]];
    for (int j=0; j<A_size[0]; j++) {
      X_ktR -= XtX(A[j], R[k]) * beta[A[j]];
    }
    
    if (std::abs(X_ktR) > lam) {
      A[A_size[0]++] = R[k]; // Add elem to A
      A_add_vec[ A_add_size_num++ ] = R[ k ];
      R[k] = R[--R_size[0]]; // Remove kth elem of R
      out=2;
    } else if (std::abs(X_ktR) >= strong_thresh) {
      S[S_size[0]++] = R[k]; // Add elem to S
      
      R[k] = R[--R_size[0]]; // Remove kth elem of R
    }
  }
  
  return out;
}

// Should be named lasso_path?
// updates beta, A, S, R
// if used in C++ should also update A_size etc. Need to check
// the first col of beta should be contain a starting value
// we should have beta.ncol() = lambda.size()
// Continues the lasso given starting values

int lasso_init(NumericMatrix beta, IntegerVector A, IntegerVector A_size,
               IntegerVector S, IntegerVector S_size, IntegerVector R, IntegerVector R_size,
               NumericVector const Xty, NumericMatrix const XtX, NumericVector const lambda,
               double const thresh, int const maxit) {
  
  for (int l=0; l<(lambda.size()-1); l++) {
    NumericMatrix::Column betacol = beta( _, l);
    do {
      beta_active(betacol, Xty, XtX, A, A_size[0], thresh, maxit, lambda[l]);
    } while (KKT_check(A, A_size, S, S_size, R, R_size, betacol, Xty, XtX,
                       lambda[l], lambda[l+1]) > 0);
    // Copy to next l for warm start
    for (int k=0; k<A_size[0]; k++) {
      beta(A[k], l+1) = betacol[A[k]];
    }
  }
  
  {
    // Brackets not strictly necessary
    int l=lambda.size()-1;
    NumericMatrix::Column betacol = beta( _, l);
    do {
      beta_active(betacol, Xty, XtX, A, A_size[0], thresh, maxit, lambda[l]);
    } while (KKT_check(A, A_size, S, S_size, R, R_size, betacol, Xty, XtX,
                       lambda[l], lambda[l]) > 0);
  }
  // no l+1 here as it will be beyond the length of lambda
  return 0;
}



// [[Rcpp::export]]
arma::sp_mat lasso_grid_sparse(arma::sp_mat beta_grid, 
                               IntegerMatrix A_add, IntegerVector A_add_size,
                               IntegerVector var_order, IntegerVector grid_seq, // rest of arguments are as usual
                               IntegerVector A, IntegerVector A_size,
                               IntegerVector S, IntegerVector S_size, IntegerVector R, IntegerVector R_size,
                               NumericVector const Xty, NumericMatrix const XtX, NumericVector const lambda,
                               double const thresh, int const maxit, double null_dev, double lambda_0, bool early_stopping )
{
  // This function would be the outer wrapper that then calls lasso_init to perform the inner regressions
  // from suitable points.
  // In order to do this we maintain record of `active sets' throughout all of the models
  IntegerVector copy_point(grid_seq.size());
  IntegerVector early_stop_point(grid_seq.size());
  int p = Xty.size();
  bool stopped_early = false;
  double eps_hat;
  int R_counter;
  std::reverse(grid_seq.begin(), grid_seq.end()); // Here we're fitting the biggest models first, other way around to lasso_var
  int prev_grid_counter, grid_counter = 0;
  int grid_size = grid_seq.size();
  int lambda_size = lambda.size();
  int start_path;
  bool new_path_point;
  double ll;
  A_size[ 0 ] = 0;
  for ( int m = 0; m < grid_size; ++m ) {
    stopped_early = false;
    A_size[ 0 ] = 0; // Build up the active set as we go along from lambda_max, for each point in the grid 
    prev_grid_counter = grid_counter;
    grid_counter = grid_seq[ m ];
    start_path = m * lambda_size;
    new_path_point = false;
    if ( m > 0 ) {
      
      checkUserInterrupt();
      while ( new_path_point == false && start_path < (m + 1) * lambda_size ) {
        for ( int k_a = 0; k_a < A_add_size[ start_path - lambda_size ]; ++k_a ) {
          for ( int k_e = prev_grid_counter; k_e < grid_counter; ++k_e ) {
            if ( A_add( k_a, start_path - lambda_size ) == var_order[ k_e ] ) {
              new_path_point = true;
            } 
          }
          if ( new_path_point == false ) {
            A_add( k_a, start_path ) = A_add( k_a, start_path - lambda_size );
            ++A_add_size[ start_path ];
          }
        }
        
        
        if ( new_path_point == false ) { 
          // Here, copy over all estimates from start_path - lambda_size, A, A_size (no need to update A_add and A_add_size)
          for ( int k_a = 0; k_a < A_add_size[ start_path - lambda_size ]; ++k_a ) {
            A[ A_size[ 0 ]++ ] = A_add( k_a, start_path - lambda_size );
          }
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            beta_grid( A[ k ], start_path ) = beta_grid( A[ k ], start_path - lambda_size );
          }
          ++start_path;
          
        } else {
          copy_point[ m ] = start_path + 1 - (m * lambda_size);
          // Rcout << "Finished copying over at path point " << start_path << std::endl;
          A_add_size[ start_path ] = 0;
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            beta_grid( A[ k ], start_path ) = beta_grid( A[ k ], start_path - 1 );
          }
        }
      }
    }
    
    
    if ( start_path < (m + 1) * lambda_size ) { // Here we call to MODIFIED lasso_init for the rest of the path
      // Now that we're going to do some computation, start treating the current estimate as a vector with entries
      // from the sparse matrix
      
      start_path = std::max(m * lambda_size, start_path - 1);
      IntegerVector range_compute = Range(start_path, (m + 1) * lambda_size - 1); // what do we do with this?
      
      S_size[ 0 ] = 0;
      R_size[ 0 ] = p - grid_counter - A_size[ 0 ];
      R_counter = 0;
      // Need to make R contain all the correct variables at the start
      // This is a source of bugs. It's where we remove all the variables that are in A from R
      for ( int k = 0; k < p - grid_counter; ++k ) { 
        R[ R_counter++ ] = var_order[ p - 1 - k ];
        for ( int kk = 0; kk < A_size[ 0 ]; ++kk ) {
          if ( A[ kk ] == R[ R_counter - 1 ] ) {
            --R_counter;
            break; //Was missing before. Need this so that we remove at most one thing from R, not more.
          } 
        }
      }
      
      // Make this before the for loop, as we can keep this vector updated throughout
      NumericVector betacol(p);
      for ( int k = 0; k < A_size[ 0 ]; ++k ) {
        betacol[ A[ k ] ] = beta_grid( A[ k ], start_path );
      }
      
      
      // Now copying the code in from lasso_init
      for ( int l = start_path; l < ((m + 1) * lambda_size - 1); ++l ) {
        checkUserInterrupt();
        ll = l - m * lambda_size;
        // This is the line that needs redoing. Store is as a simple dense NumericVector
        
        // NumericMatrix::Column betacol = beta_grid(  _, l );
        if ( stopped_early == false ) {
          IntegerMatrix::Column A_add_vec = A_add( _, l );
          do {
            beta_active(betacol, Xty, XtX, A, A_size[0], thresh, maxit, lambda[ ll ]);
          } while (KKT_check_m(A_add_vec, A_add_size[ l ], 
                               A, A_size, S, S_size, R, R_size, betacol, Xty, XtX,
                               lambda[ ll ], lambda[ ll + 1 ]) > 0);
          
          // Input the current solution vector into the sparse matrix grid. A is only increasing so this overwrites
          // all of the non-zero entries from the warm start
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            beta_grid( A[ k ], l ) = betacol[ A[ k ] ];
          }
          
          // // Copy to next l for warm start. This isn't actually necessary now that we carry over betacol from before!
          // // Dare we remove it?
          // for (int k=0; k<A_size[0]; k++) {
          //   beta_grid(A[k], l+1) = betacol[A[k]];
          // }
          
          for ( int k_a = 0; k_a < A_add_size[ l ]; ++k_a ) {
            A_add( k_a, l ) = A_add_vec[ k_a ];
          }
          
          
          
          
          // End of for loop for computing answer
          // Compute square root objective value here and see if early termination good
          eps_hat = null_dev;
          for ( int k_a = 0; k_a < A_size[ 0 ]; ++k_a ) {
            eps_hat -= 2 * betacol[ A[ k_a ] ] * Xty[ A[ k_a ] ];
            for ( int k_a2 = 0; k_a2 < A_size[ 0 ]; ++k_a2 ) {
              eps_hat += betacol[ A[ k_a ] ] * XtX( k_a, k_a2 ) * betacol[ A[ k_a2 ] ];
            }
          }
          
          // eps_hat should now be the value of 1/n || Y - X \beta ||_2^2
          if ( eps_hat * pow(lambda_0, 2) - pow(lambda[ ll ], 2) > 0 ) {
            stopped_early = true; // 9Oct
            early_stop_point[ m ] = ll + 1;
            // Rcout << "stopping early criterion met; m = " << m << ", and lambda = " << ll << std::endl;
            if ( early_stopping == false ) stopped_early = false;
            // Rcout << "Stopped early at gridpoint " << ll << " at set " << m << std::endl;
          }
        } else {
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            // beta_grid( A[ k ], l ) = beta_grid( A[ k ], l - 1 );
            beta_grid( A[ k ], l ) = betacol[ A[ k ] ];
            // Rcout << "Stopped early at girdpoint " << ll << " at set " << m << std::endl;
          }
        } 
      }
      
      {
        // Brackets not strictly necessary
        checkUserInterrupt();
        int l = (m + 1) * lambda_size - 1;
        ll = l - m * lambda_size;
        if ( stopped_early == false ) {
          // NumericMatrix::Column betacol = beta_grid( _, l); // Again don't need this in the sparse matrix setup
          IntegerMatrix::Column A_add_vec = A_add( _, l );
          do {
            beta_active(betacol, Xty, XtX, A, A_size[0], thresh, maxit, lambda[ ll ]);
          } while (KKT_check_m(A_add_vec, A_add_size[ l ], 
                               A, A_size, S, S_size, R, R_size, betacol, Xty, XtX,
                               lambda[ ll ], lambda[ ll ]) > 0);
          
          for ( int k_a = 0; k_a < A_add_size[ l ]; ++k_a ) {
            A_add( k_a, l ) = A_add_vec[ k_a ];
          }
          
          // Now copy in final beta column into the grid
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            beta_grid( A[ k ], l ) = betacol[ A[ k ] ];
          }
        } else {
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            
            // beta_grid( A[ k ], l ) = beta_grid( A[ k ], l - 1 );
            beta_grid( A[ k ], l ) = betacol[ A[ k ] ];
            
            // beta_grid( A[ k ], l ) = beta_grid( A[ k ], l - 1 );
            
          }
          
        }
        // no l+1 here as it will be beyond the length of lambda
        
      }
    }
  }
  // Can also get it to return copy_point and early_stop_point for diagnostics of how early the paths are being copied / terminated under early stopping.
  return beta_grid;
}


// [[Rcpp::export]]
int lasso_grid(NumericMatrix beta_grid, 
               IntegerMatrix A_add, IntegerVector A_add_size,
               IntegerVector var_order, IntegerVector grid_seq, // rest of arguments are as usual
               IntegerVector A, IntegerVector A_size,
               IntegerVector S, IntegerVector S_size, IntegerVector R, IntegerVector R_size,
               NumericVector const Xty, NumericMatrix const XtX, NumericVector const lambda,
               double const thresh, int const maxit)
{
  // This function would be the outer wrapper that then calls lasso_init to perform the inner regressions
  // from suitable points.
  // In order to do this we maintain record of `active sets' throughout all of the models
  int p = Xty.size();
  int R_counter;
  std::reverse(grid_seq.begin(), grid_seq.end()); // Here we're fitting the biggest models first, other way around to lasso_var
  int prev_grid_counter, grid_counter = 0;
  int grid_size = grid_seq.size();
  int lambda_size = lambda.size();
  int start_path;
  bool new_path_point;
  double ll;
  A_size[ 0 ] = 0;
  for ( int m = 0; m < grid_size; ++m ) {
    A_size[ 0 ] = 0; // Build up the active set as we go along from lambda_max, for each point in the grid 
    prev_grid_counter = grid_counter;
    grid_counter = grid_seq[ m ];
    start_path = m * lambda_size;
    new_path_point = false;
    
    if ( m > 0 ) {
      while ( new_path_point == false && start_path < (m + 1) * lambda_size ) {
        for ( int k_a = 0; k_a < A_add_size[ start_path - lambda_size ]; ++k_a ) {
          for ( int k_e = prev_grid_counter; k_e < grid_counter; ++k_e ) {
            if ( A_add( k_a, start_path - lambda_size ) == var_order[ k_e ] ) {
              new_path_point = true;
            } 
          }
          if ( new_path_point == false ) {
            A_add( k_a, start_path ) = A_add( k_a, start_path - lambda_size );
            ++A_add_size[ start_path ];
          }
        }
        
        
        if ( new_path_point == false ) { 
          // Here, copy over all estimates from start_path - lambda_size, A, A_size (no need to update A_add and A_add_size)
          for ( int k_a = 0; k_a < A_add_size[ start_path - lambda_size ]; ++k_a ) {
            A[ A_size[ 0 ]++ ] = A_add( k_a, start_path - lambda_size );
          }
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            beta_grid( A[ k ], start_path ) = beta_grid( A[ k ], start_path - lambda_size );
          }
          ++start_path;
          
        } else {
          A_add_size[ start_path ] = 0;
          for ( int k = 0; k < A_size[ 0 ]; ++k ) {
            beta_grid( A[ k ], start_path ) = beta_grid( A[ k ], start_path - 1 );
          }
        }
      }
    }
    
    
    if ( start_path < (m + 1) * lambda_size ) { // Here we call to MODIFIED lasso_init for the rest of the path
      start_path = std::max(m * lambda_size, start_path - 1);
      IntegerVector range_compute = Range(start_path, (m + 1) * lambda_size - 1);
      
      S_size[ 0 ] = 0;
      R_size[ 0 ] = p - grid_counter - A_size[ 0 ];
      R_counter = 0;
      // Need to make R contain all the correct variables at the start
      for ( int k = 0; k < p - grid_counter; ++k ) { 
        R[ R_counter++ ] = var_order[ p - 1 - k ];
        for ( int kk = 0; kk < A_size[ 0 ]; ++kk ) {
          if ( A[ kk ] == R[ R_counter - 1 ] ) {
            --R_counter;
          } 
        }
      }
      
      
      // Now copying the code in from lasso_init
      for ( int l = start_path; l < ((m + 1) * lambda_size - 1); ++l ) {
        ll = l - m * lambda_size;
        NumericMatrix::Column betacol = beta_grid(  _, l );
        IntegerMatrix::Column A_add_vec = A_add( _, l );
        do {
          beta_active(betacol, Xty, XtX, A, A_size[0], thresh, maxit, lambda[ ll ]);
        } while (KKT_check_m(A_add_vec, A_add_size[ l ], 
                             A, A_size, S, S_size, R, R_size, betacol, Xty, XtX,
                             lambda[ ll ], lambda[ ll + 1 ]) > 0);
        // Copy to next l for warm start
        
        
        for (int k=0; k<A_size[0]; k++) {
          
          beta_grid(A[k], l+1) = betacol[A[k]];
        }
        
        for ( int k_a = 0; k_a < A_add_size[ l ]; ++k_a ) {
          A_add( k_a, l ) = A_add_vec[ k_a ];
        }
        
      }
      
      {
        // Brackets not strictly necessary
        int l = (m + 1) * lambda_size - 1;
        ll = l - m * lambda_size;
        NumericMatrix::Column betacol = beta_grid( _, l);
        IntegerMatrix::Column A_add_vec = A_add( _, l );
        ////21SEP Rcout << "l_g 12" << std::endl;
        do {
          beta_active(betacol, Xty, XtX, A, A_size[0], thresh, maxit, lambda[ ll ]);
        } while (KKT_check_m(A_add_vec, A_add_size[ l ], 
                             A, A_size, S, S_size, R, R_size, betacol, Xty, XtX,
                             lambda[ ll ], lambda[ ll ]) > 0);
        
        for ( int k_a = 0; k_a < A_add_size[ l ]; ++k_a ) {
          A_add( k_a, l ) = A_add_vec[ k_a ];
        }
        
      }
      // no l+1 here as it will be beyond the length of lambda
      
    }
  }
  
  return 0;
}

