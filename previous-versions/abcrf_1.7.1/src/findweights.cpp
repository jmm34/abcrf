#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix findweights(NumericMatrix origNodes, IntegerMatrix inbag, NumericMatrix nodes, int nobs, int nnew, int ntree){
  NumericMatrix result(nobs,nnew);
  IntegerVector counti(nobs);
  double meancount;
  for(int k=0; k<ntree; k++){
   for(int i=0; i<nnew; i++){
     meancount = 0;
     for(int j=0; j<nobs; j++){
       if( ( origNodes(j,k) == nodes(i,k) ) && (inbag(j,k) != 0) ){
         counti[j] = inbag(j,k);
         meancount = meancount + inbag(j,k);
       } else{
         counti[j] = 0;
       }
     }
     if( meancount >=1 ){
       for(int j=0; j<nobs; j++){
         result(j,i) = result(j,i) + counti[j]/ meancount;
       }
     }
   } 
  }
  return(result);
}
