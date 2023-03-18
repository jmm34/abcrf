// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// findweights
NumericMatrix findweights(NumericMatrix origNodes, IntegerMatrix inbag, NumericMatrix nodes, int nobs, int nnew, int ntree);
RcppExport SEXP abcrf_findweights(SEXP origNodesSEXP, SEXP inbagSEXP, SEXP nodesSEXP, SEXP nobsSEXP, SEXP nnewSEXP, SEXP ntreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type origNodes(origNodesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type inbag(inbagSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type nnew(nnewSEXP);
    Rcpp::traits::input_parameter< int >::type ntree(ntreeSEXP);
    rcpp_result_gen = Rcpp::wrap(findweights(origNodes, inbag, nodes, nobs, nnew, ntree));
    return rcpp_result_gen;
END_RCPP
}
// oobErrors
NumericVector oobErrors(IntegerVector sequo, int ntrain, IntegerVector mod, int ntree, IntegerVector modindex, IntegerMatrix inbag, IntegerMatrix mimi);
RcppExport SEXP abcrf_oobErrors(SEXP sequoSEXP, SEXP ntrainSEXP, SEXP modSEXP, SEXP ntreeSEXP, SEXP modindexSEXP, SEXP inbagSEXP, SEXP mimiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type sequo(sequoSEXP);
    Rcpp::traits::input_parameter< int >::type ntrain(ntrainSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mod(modSEXP);
    Rcpp::traits::input_parameter< int >::type ntree(ntreeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type modindex(modindexSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type inbag(inbagSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type mimi(mimiSEXP);
    rcpp_result_gen = Rcpp::wrap(oobErrors(sequo, ntrain, mod, ntree, modindex, inbag, mimi));
    return rcpp_result_gen;
END_RCPP
}
// oobErrorsReg
NumericVector oobErrorsReg(IntegerVector sequo, int ntrain, int ntree, NumericVector resp, IntegerMatrix inbag, NumericMatrix pred);
RcppExport SEXP abcrf_oobErrorsReg(SEXP sequoSEXP, SEXP ntrainSEXP, SEXP ntreeSEXP, SEXP respSEXP, SEXP inbagSEXP, SEXP predSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type sequo(sequoSEXP);
    Rcpp::traits::input_parameter< int >::type ntrain(ntrainSEXP);
    Rcpp::traits::input_parameter< int >::type ntree(ntreeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type resp(respSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type inbag(inbagSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pred(predSEXP);
    rcpp_result_gen = Rcpp::wrap(oobErrorsReg(sequo, ntrain, ntree, resp, inbag, pred));
    return rcpp_result_gen;
END_RCPP
}
// predictOob_cpp
NumericMatrix predictOob_cpp(NumericMatrix origNodes, NumericMatrix inbag, int nobs, int ntree);
RcppExport SEXP abcrf_predictOob_cpp(SEXP origNodesSEXP, SEXP inbagSEXP, SEXP nobsSEXP, SEXP ntreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type origNodes(origNodesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type inbag(inbagSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type ntree(ntreeSEXP);
    rcpp_result_gen = Rcpp::wrap(predictOob_cpp(origNodes, inbag, nobs, ntree));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"abcrf_findweights", (DL_FUNC) &abcrf_findweights, 6},
    {"abcrf_oobErrors", (DL_FUNC) &abcrf_oobErrors, 7},
    {"abcrf_oobErrorsReg", (DL_FUNC) &abcrf_oobErrorsReg, 6},
    {"abcrf_predictOob_cpp", (DL_FUNC) &abcrf_predictOob_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_abcrf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}