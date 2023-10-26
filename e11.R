library("iterators")
library("parallel")
library("foreach")
library("doParallel")
path="/nethome/jxz280/neuron/power/erdos1/"
source(paste(path,"erdos.R",sep=""))
set.seed(1)
######### data generation ####################
registerDoParallel(cores=8)
ht1=foreach(m=1:50) %dopar%{
  cat("simulation",m,"\n")
  onesim(84,20,netstr=net.use1)
}
save(ht1,file=paste(path,"ht1.RData",sep=""))

######### estimation ##########################
registerDoParallel(cores=8)
oneest1=foreach(m=1:50) %dopar%{
  cat("simulation",m,"\n")
  h0=ht1[[m]]
  G1=G.est(h0,20,8)
  b0.all=sapply(1:84,function(i){
    cat("neuron",i,"\n")
    b.est(h0,19.999,i,n1=8)
  })
  list(G0=G1,b.all=b0.all)
}
save(oneest1,file=paste(path,"oneest1.RData",sep=""))

######### tunning ##############################
registerDoParallel(cores=8)
a.all.total1=foreach(m=1:50) %dopar%{
  cat("simulation",m,"\n")
  G0=oneest1[[m]]$G0/20
  a.all=lapply(1:84,function(i){
    b0=oneest1[[m]]$b.all[,i]/20
    
    lambda0=exp(seq(0,4,length=100))
    a.total= sapply(1:length(lambda0),function(j){
      lambda=lambda0[j]
      grplas0(G0,b0,lambda,M=84,n1=8,n0=20)
    })
    
    #log(T)log(p)^2
    loss=sapply(1:100,function(j){
      a0=a.total[,j]
      nonzero=sum(sapply(1:84,function(k) sum(a0[(8*k+13):(8*k+20)]==0)==0))
      2*(-2*t(a0)%*%b0+t(a0)%*%G0%*%a0)*(20)/sum(b0[1:20])+nonzero*log(20)*(log(84))^2
    } )
    id0=which.min(loss)
    a0=a.total[,id0]
    
    list(a0=a0,a.total=a.total,eta=lambda0[id0])
  })
  a.all
}
save(a.all.total1,file=paste(path,"a.all.total1.RData",sep=""))

rate1 <- sapply(1:50,function(m){
  a.all=a.all.total1[[m]]
  a=sapply(1:84,function(i)a.all[[i]]$a0)
  dir=sapply(1:84,function(i)sapply(1:84,function(j){
    a[8*j+13,i]
  }))
  dir[dir>0]=1
  dir[dir<0]=-1
  rate=econ3(dir,netstr=net.use1)
  c(rate$tfnr,rate$fpr,rate$f1)
})
save(rate1,file=paste(path,"rate1.RData",sep=""))