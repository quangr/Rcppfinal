#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// [[Rcpp::export]]
mat func2(mat X,int m){
  int c_m2 = m*(m-1)/2;
  mat combinematrix(c_m2,2);
  int i = 0;
  for(int j = 1;j <= m;j++){
    for(int k = j+1;k <= m;k++){
      combinematrix(i,0) = j;
      combinematrix(i,1) = k;
      i++;
    }
  }
  
  mat amatrix(c_m2,m, fill::zeros);
  for(int l = 0;l < c_m2;l++){
    amatrix(l,combinematrix(l,0)-1) = 1;
    amatrix(l,combinematrix(l,1)-1) = -1;
  }

  mat resultmatrix(X.n_cols*c_m2,X.n_rows*m);
  mat temp(c_m2,m);
  for(int i = 0;i<X.n_rows;i++){
    for(int j = 0;j<X.n_cols;j++){
      temp = amatrix * X(i,j);
      for(int k = 0;k<c_m2;k++){
        for(int l = 0;l<m;l++){
          resultmatrix(j*c_m2+k,i*m+l) = temp(k,l);
        }
      }
    }
  }
  
  return resultmatrix;
}