library(ROI)
library(purrr)

# constraint matrix 1, input:X(n*r), return a constraint matrix A((p(p-1)r/2)*(np))
PSOD.Finite.qnorm.1.constraint1 <- function(X, m) {
  combo <- map(1:(m - 1), ~ cross2(c(.x), (.x + 1):m)) %>% flatten()
  print(combo)
  commatrix <- map(combo, function(v) {
    t <- rep(0, m)
    t[v[[1]]] <- 1
    t[v[[2]]] <- (-1)
    t%>%matrix(1,m)
  }) %>% reduce(rbind)
  sintraitmatrix <- function(x) {
    map(x, ~ .x * commatrix) %>% reduce(cbind)
  }
  map(as.data.frame(X), sintraitmatrix) %>% reduce(rbind)%>%unname()
}


# constraint matrix 1, input:X(n*r), return a constraint matrix A((2^r)*(np))
PSOD.Finite.qnorm.inf.constraint1 <- function(X, m) {
  combo <- map(1:(m - 1), ~ cross2(c(.x), (.x + 1):m)) %>% flatten()
  print(combo)
  commatrix <- map(combo, function(v) {
    t <- rep(0, m)
    t[v[[1]]] <- 1
    t[v[[2]]] <- (-1)
    t%>%matrix(1,m)
  }) %>% reduce(rbind)
  sintraitmatrix <- function(x) {
    map(x, ~ .x * commatrix) %>% reduce(cbind)
  }
  if(ncol(X)==1){
    crossones=list(list(-1,1))
  }else{
    crossones=reduce(rep(1,ncol(X)-1),~map(c(-1,1),function(x){cross2(c(x),.x)%>%map(flatten)})%>%flatten(),.init =c(-1,1))
  }
  matrixs=map(as.data.frame(X), sintraitmatrix)
  map(crossones,~.x%>%reduce2(matrixs,function(x,y,z){x+y*z},.init = 0))%>%reduce(rbind)
}

PSOD.common.constraint1 <- function(n, m) {
  map(1:n,~c(rep(0,(.x-1)*m),rep(1,m),rep(0,(n-.x)*m))%>%reduce(cbind))%>%reduce(rbind)%>%unname()
}

PSOD.common.constraint2 <- function(n, m) {
  map(1:n,~diag(m))%>%reduce(cbind)%>%unname()
}

# Finite dimensional q-space
PSOD.Finite.qnorm <- function(data, q = 1, formula, m = 2) {
  design <- model.matrix(formula, data = data)
  n=nrow(design)
  
  CCL1=PSOD.common.constraint1(n,m)
  CCL1=cbind(rep(0,nrow(CCL1)),CCL1)
  CCL2=PSOD.common.constraint2(n,m)
  CCL2=cbind(rep(0,nrow(CCL2)),CCL2)
  CCL1=L_constraint(CCL1,dir=rep("==",nrow(CCL1)),rhs=rep(1,nrow(CCL1)))
  CCL2=L_constraint(CCL2,dir=rep("==",nrow(CCL2)),rhs=rep(n/m,nrow(CCL2)))
  
  if (q == 1) {
    CL1=PSOD.Finite.qnorm.1.constraint1(design,m)
    CL1=rbind(CL1,-CL1)
    CL1=cbind(rep(1,nrow(CL1)),CL1)
    CL1=L_constraint(CL1,dir=rep(">=",nrow(CL1)),rhs=rep(0,nrow(CL1)))

        types = c("C",rep("B", m*n))
    milp <- OP(objective = L_objective(c(1, rep(0, m*n))),
               constraints = c(CL1,CCL1,CCL2),
               types =types, 
               maximum = FALSE)
    (sol <- ROI_solve(milp))
    solution(sol)
  }

  if (q == Inf) {
    CL1=PSOD.common.constraint1(design,m)
    CL1=rbind(CL1,-CL1)
    CL1=cbind(rep(1,nrow(CL1)),CL1)
    CL1=L_constraint(CL1,dir=rep(">=",nrow(CL1)),rhs=rep(0,nrow(CL1)))
    
    types = c("C",rep("B", m*n))
    milp <- OP(objective = L_objective(c(1, rep(0, m*n))),
               constraints = c(CL1,CCL1,CCL2),
               types =types, 
               maximum = FALSE)
    (sol <- ROI_solve(milp))
    solution(sol)
  }
  
  
}

# Lipschitz functions


PSOD.Lipschitz <- function(data, X, q = 1) {

}


# Reproducing kernel Hilbert space


PSOD.RKHS <- function(y, X, q = 1) {

}
