library(ROI)
library(purrr)
library(furrr)
future::plan(multisession)
source("./PSOD.R", local = T)
RunSimulation=function(){
  future_map(1:10,function(i){
    test=data.frame(x=rnorm(10),y=rnorm(10))
    test=test%>%rbind(test+data.frame(x=rnorm(10,sd = 0.05),y=rnorm(10,sd = 0.05)))
    ar1=matrix(PSOD.Finite.qnorm(test),ncol = 2,byrow = T)%*%c(1,2)%>%as.vector()
    ar2=matrix(PSOD.Finite.qnorm(test,q=Inf),ncol = 2,byrow = T)%*%c(1,2)%>%as.vector()
    ar3=rep(c(1,2),each=10)
    treat1=test$x-test$y+1/2
    treat2=test$x-test$y-1/2
    map(list(ar1,ar2,ar3),~mean(treat1[.x==1]-treat2[.x==2]))
  })%>%transpose()%>%map(~var(unlist(.x)))%>%unlist
  
}