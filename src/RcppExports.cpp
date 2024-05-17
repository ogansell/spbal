// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cppBASpts
Rcpp::List cppBASpts(int n, Rcpp::IntegerVector seeds, Rcpp::NumericVector bases, bool verbose);
RcppExport SEXP _spbal_cppBASpts(SEXP nSEXP, SEXP seedsSEXP, SEXP basesSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bases(basesSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cppBASpts(n, seeds, bases, verbose));
    return rcpp_result_gen;
END_RCPP
}
// cppRSHalton_br
Rcpp::List cppRSHalton_br(int n, Rcpp::NumericVector bases, Rcpp::NumericVector seeds, bool verbose);
RcppExport SEXP _spbal_cppRSHalton_br(SEXP nSEXP, SEXP basesSEXP, SEXP seedsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bases(basesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cppRSHalton_br(n, bases, seeds, verbose));
    return rcpp_result_gen;
END_RCPP
}
// cppBASptsIndexed
Rcpp::List cppBASptsIndexed(int n, Rcpp::IntegerVector seeds, Rcpp::NumericVector bases, Rcpp::IntegerVector boxes, bool verbose);
RcppExport SEXP _spbal_cppBASptsIndexed(SEXP nSEXP, SEXP seedsSEXP, SEXP basesSEXP, SEXP boxesSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bases(basesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type boxes(boxesSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cppBASptsIndexed(n, seeds, bases, boxes, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spbal_cppBASpts", (DL_FUNC) &_spbal_cppBASpts, 4},
    {"_spbal_cppRSHalton_br", (DL_FUNC) &_spbal_cppRSHalton_br, 4},
    {"_spbal_cppBASptsIndexed", (DL_FUNC) &_spbal_cppBASptsIndexed, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_spbal(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}