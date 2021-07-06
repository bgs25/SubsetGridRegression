// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lasso_grid_sparse
arma::sp_mat lasso_grid_sparse(arma::sp_mat beta_grid, IntegerMatrix A_add, IntegerVector A_add_size, IntegerVector var_order, IntegerVector grid_seq, IntegerVector A, IntegerVector A_size, IntegerVector S, IntegerVector S_size, IntegerVector R, IntegerVector R_size, NumericVector const Xty, NumericMatrix const XtX, NumericVector const lambda, double const thresh, int const maxit, double null_dev, double lambda_0, bool early_stopping);
RcppExport SEXP _SubsetGridRegression_lasso_grid_sparse(SEXP beta_gridSEXP, SEXP A_addSEXP, SEXP A_add_sizeSEXP, SEXP var_orderSEXP, SEXP grid_seqSEXP, SEXP ASEXP, SEXP A_sizeSEXP, SEXP SSEXP, SEXP S_sizeSEXP, SEXP RSEXP, SEXP R_sizeSEXP, SEXP XtySEXP, SEXP XtXSEXP, SEXP lambdaSEXP, SEXP threshSEXP, SEXP maxitSEXP, SEXP null_devSEXP, SEXP lambda_0SEXP, SEXP early_stoppingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type beta_grid(beta_gridSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type A_add(A_addSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A_add_size(A_add_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type var_order(var_orderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type grid_seq(grid_seqSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A_size(A_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type S_size(S_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type R(RSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type R_size(R_sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< NumericMatrix const >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double const >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int const >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type null_dev(null_devSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_0(lambda_0SEXP);
    Rcpp::traits::input_parameter< bool >::type early_stopping(early_stoppingSEXP);
    rcpp_result_gen = Rcpp::wrap(lasso_grid_sparse(beta_grid, A_add, A_add_size, var_order, grid_seq, A, A_size, S, S_size, R, R_size, Xty, XtX, lambda, thresh, maxit, null_dev, lambda_0, early_stopping));
    return rcpp_result_gen;
END_RCPP
}
// lasso_grid
int lasso_grid(NumericMatrix beta_grid, IntegerMatrix A_add, IntegerVector A_add_size, IntegerVector var_order, IntegerVector grid_seq, IntegerVector A, IntegerVector A_size, IntegerVector S, IntegerVector S_size, IntegerVector R, IntegerVector R_size, NumericVector const Xty, NumericMatrix const XtX, NumericVector const lambda, double const thresh, int const maxit);
RcppExport SEXP _SubsetGridRegression_lasso_grid(SEXP beta_gridSEXP, SEXP A_addSEXP, SEXP A_add_sizeSEXP, SEXP var_orderSEXP, SEXP grid_seqSEXP, SEXP ASEXP, SEXP A_sizeSEXP, SEXP SSEXP, SEXP S_sizeSEXP, SEXP RSEXP, SEXP R_sizeSEXP, SEXP XtySEXP, SEXP XtXSEXP, SEXP lambdaSEXP, SEXP threshSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta_grid(beta_gridSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type A_add(A_addSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A_add_size(A_add_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type var_order(var_orderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type grid_seq(grid_seqSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A_size(A_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type S_size(S_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type R(RSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type R_size(R_sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< NumericMatrix const >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double const >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int const >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(lasso_grid(beta_grid, A_add, A_add_size, var_order, grid_seq, A, A_size, S, S_size, R, R_size, Xty, XtX, lambda, thresh, maxit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SubsetGridRegression_lasso_grid_sparse", (DL_FUNC) &_SubsetGridRegression_lasso_grid_sparse, 19},
    {"_SubsetGridRegression_lasso_grid", (DL_FUNC) &_SubsetGridRegression_lasso_grid, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_SubsetGridRegression(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}