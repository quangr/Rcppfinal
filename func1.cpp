#include <iostream>
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;

vector<vector<int>> ordered_pair(int m)
{
  vector<vector<int>> res;
  for (int i=1; i <= (m-1); i++)
  {
    for (int j = i+1; j <= m; j++)
    {
      vector<int> tmp = {i,j};
      res.push_back(tmp);
    }

  }
  return res;

}

vector<int> matrix_realligned1(mat A)
{
  vector<int> res;
  int n = A.n_cols;
  for (int i=0; i <= (n-1); i++)
  {
    for (int j=0; j <= (n-1); j++)
    {
      res.push_back(A.at(i,j));
    }
  }
  return res;
}

vec matrix_realligned(mat A)
{
  int n = A.n_cols;
  vec res(n*n);
  for (int i=0; i <= (n-1); i++)
  {
    for (int j=0; j <= (n-1); j++)
    {
      res.at(i*n + j) = A.at(i,j);
    }
  }
  return res;
}

// [[Rcpp::export]]
mat PSOD_Lipschitz_constraint(int n, int m)
{
  vector<vector<int>> pair = ordered_pair(m);
  mat M;
  for (int k=1; k<= pair.size(); k++)
  {
    int a = pair.at(k-1).at(0);
    int b = pair.at(k-1).at(1);
    vec l(m);
    l.at(a-1) = -1;
    l.at(b-1) = 1;
    mat Mk(n, m*n + n*n*m*(m-1)/2, fill::zeros);
    for (int i=1; i<=n; i++)
    {
      vec vi(m*n, fill::zeros);
      for (int t=m*(i-1); t <= m*(i-1)+(m-1); t++)
      {
        vi.at(t) = l.at(t-m*(i-1));
      }
      vec ei(n);
      ei.at(i-1) = 1;
      vec vector1(n, fill::ones);
      mat Ai = vector1 * ei.t() - ei * vector1.t();
      vec lazhi = matrix_realligned(Ai);
      vec ui(n*n*m*(m-1)/2, fill::zeros);
      for (int t=n*n*(k-1); t <= n*n*(k-1) + (n*n-1); t++)
      {
        ui.at(t) = lazhi.at(t-n*n*(k-1));
      }
      for (int j=0; j <= m*n-1; j++)
      {
        Mk.at(i-1,j) = vi.at(j);
      }
      
      for (int j = 0; j<= n*n*m*(m-1)/2-1; j++)
      {
        Mk.at(i-1, j+m*n) = ui.at(j);
      }
    }
    
    if (k == 1){M = Mk;}
    else{ M = join_cols(M, Mk);}
  }
  return M;
}
