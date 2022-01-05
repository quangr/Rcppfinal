#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// [[Rcpp::export]]
mat fun3(int k,mat x){
  int n=k*(k-1)/2;
  mat order(n,2);
  int kk=0;
  for(int i=1;i<k+1;i+=1){
    for(int j=i+1;j<k+1;j+=1){
      order(kk,0)=i;
      order(kk,1)=j;
      kk+=1;
    }
  }
  mat A(n,k,fill::zeros);
  for(int i=0;i<n;i+=1){
    A(i,order(i,0)-1)=1;
    A(i,order(i,1)-1)=-1;
  }
  int xn=x.n_cols;
  cube c(n,n*n,xn);
  for(int i=0;i<xn;i+=1){
    mat v=x(0,i)*a;
    for(int j=1;j<n;j+=1){
      v=join_rows(v,x(j,i)*a);
    }
    c.slice(i)=v;
  }
  mat re(n,n*n,fill::zeros);
  for(int i=0;i<xn;i+=1){
    re-=c.slice(i);
  }
  for(int i=1;i<pow(2,xn);i++){
    mat temp(n,n*n,fill::zeros);
    int num=i;
    for(int j=0;j<xn;j++){
      if (num%2==1){temp+=c.slice(xn-j-1);}
      else {temp-=c.slice(xn-j-1);}
      num/=2;
    }
    re=join_cols(re,temp);
  }
  return re;


}


