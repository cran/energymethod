// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// energy_method
List energy_method(arma::cube sample_1, arma::cube sample_2, int num_bootstrap_reps, int seed, std::string type);
RcppExport SEXP _energymethod_energy_method(SEXP sample_1SEXP, SEXP sample_2SEXP, SEXP num_bootstrap_repsSEXP, SEXP seedSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type sample_1(sample_1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sample_2(sample_2SEXP);
    Rcpp::traits::input_parameter< int >::type num_bootstrap_reps(num_bootstrap_repsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(energy_method(sample_1, sample_2, num_bootstrap_reps, seed, type));
    return rcpp_result_gen;
END_RCPP
}
// energy_method_complex
List energy_method_complex(arma::cx_cube sample_1, arma::cx_cube sample_2, int num_bootstrap_reps, int seed, std::string type);
RcppExport SEXP _energymethod_energy_method_complex(SEXP sample_1SEXP, SEXP sample_2SEXP, SEXP num_bootstrap_repsSEXP, SEXP seedSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_cube >::type sample_1(sample_1SEXP);
    Rcpp::traits::input_parameter< arma::cx_cube >::type sample_2(sample_2SEXP);
    Rcpp::traits::input_parameter< int >::type num_bootstrap_reps(num_bootstrap_repsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(energy_method_complex(sample_1, sample_2, num_bootstrap_reps, seed, type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_energymethod_energy_method", (DL_FUNC) &_energymethod_energy_method, 5},
    {"_energymethod_energy_method_complex", (DL_FUNC) &_energymethod_energy_method_complex, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_energymethod(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
