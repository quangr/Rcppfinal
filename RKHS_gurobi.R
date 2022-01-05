
library(gurobi)
library(Matrix)



PSOD.RKHS <-  function(data,m = 2,rkfunc=rk) {
  grammat=gram(data,function(s,t){exp(-norm(s-t, type="2")) })
  
  n=nrow(grammat)
  
  PSOD.common.constraint1 <- function(n, m) {
    map(1:n,~c(rep(0,(.x-1)*m),rep(1,m),rep(0,(n-.x)*m))%>%reduce(cbind))%>%reduce(rbind)%>%unname()
  }
  
  PSOD.common.constraint2 <- function(n, m) {
    map(1:n,~diag(m))%>%reduce(cbind)%>%unname()
  }
  
  CCL1=PSOD.common.constraint1(m,n)
  CCL1=cbind(rep(0,nrow(CCL1)),CCL1)
  CCL2=PSOD.common.constraint2(m,n)
  CCL2=cbind(rep(0,nrow(CCL2)),CCL2)
  
  
  model <- list()
  
  model$A          <- rbind(CCL1,CCL2)
  model$modelsense <- 'min'
  model$obj        <- c(1, rep(0, m*n))
  model$rhs        <- c(rep(n/m,nrow(CCL1)),rep(1,nrow(CCL2)))
  model$sense      <- c('=')
  model$vtype      <- c("C",rep("B", m*n))
  combo <- map(1:(m - 1), ~ cross2(c(.x), (.x + 1):m)) %>% flatten()
  # First quadratic constraint: x^2 + y^2 - z^2 <= 0
  
  model$quadcon <- map(combo,function(v){
    qc <- list()
    qc$Qc <- bdiag(append(0,map(1:m,function(i){if(i==v[1]){return(grammat)}
      if(i==v[2]){return(-grammat)}
      diag(0,n)
    })))
    qc$q<-c(-1,rep(0, m*n))
    qc$rhs <- 0.0
    qc
  })
  
  
  
  result <- gurobi(model)
  
#  print(result$objval)
 # print(result$x)
  return(matrix(result$x[-1],m,byrow = T)%>%t()%*%1:m%>%as.vector())
}


rk <- function(s, t) {
  s%*%t
} 

gram <- function(X, rkfunc = rk) {
  apply(X, 1, function(Row) 
    apply(X, 1, function(tRow) rkfunc(Row, tRow))
  )
}

data=pairdata(20,m)

