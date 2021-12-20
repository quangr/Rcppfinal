library(ROI)
library(purrr)
library(furrr)
library(ggplot2)
future::plan(multisession)
source("./PSOD.R", local = T)
pairdata=function(n,m=2){
  test=data.frame(x=rnorm(n),y=rnorm(n))
  map(1:(m-1),~test+data.frame(x=rnorm(n,sd = 0.05),y=rnorm(n,sd = 0.05)))%>%reduce(rbind,.init = test)
}

RunSimulation=function(){
  future_map(1:10,function(i){
    test=data.frame(x=rnorm(10),y=rnorm(10))
    test=test%>%rbind(test+data.frame(x=rnorm(10,sd = 0.05),y=rnorm(10,sd = 0.05)))
    ar1=matrix(PSOD.Finite.qnorm(test),ncol = 2,byrow = T)%*%1:2%>%as.vector()
    ar2=matrix(PSOD.Finite.qnorm(test,q=Inf),ncol = 2,byrow = T)%*%1:2%>%as.vector()
    ar3=rep(c(1,2),each=10)
    treat1=test$x-test$y+1/2
    treat2=test$x-test$y-1/2
    map(list(ar1,ar2,ar3),~mean(treat1[.x==1]-treat2[.x==2]))
  })%>%transpose()%>%map(~var(unlist(.x)))%>%unlist
  
}


PlotResult=function(m){
  data=pairdata(30,m)
  ar1=matrix(PSOD.Lipschitz(data,m=m),ncol = m,byrow = T)%*%1:m%>%as.vector()
  ggplot(data,aes(x=x,y=y,shape=as.factor(ar1)))+geom_point()
}