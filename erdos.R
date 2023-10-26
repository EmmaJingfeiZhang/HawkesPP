#########################################
############ generation #################
#########################################
load(paste(path,"set_erdos/net.use1.RData",sep=""))
load(paste(path,"set_erdos/level1.RData",sep=""))
load(paste(path,"set_erdos/level2.RData",sep=""))
load(paste(path,"set_erdos/level3.RData",sep=""))
load(paste(path,"set_erdos/level4.RData",sep=""))
load(paste(path,"set_erdos/level5.RData",sep=""))
load(paste(path,"set_erdos/level6.RData",sep=""))
load(paste(path,"set_erdos/alpha.RData",sep=""))
gen0<-function(base,T,lambda.m){#inhomogeneous poisson
  points=vector()
  s=n=0
  while(s<T){
    #lambda.m=lambda.up(base,interac,points,s,T,T1)
    u=runif(1)
    w=-log(u)/lambda.m
    s=s+w
    D=runif(1)
    lambda.m0=base(s)
    if(D*lambda.m<=lambda.m0){
      n=n+1
      points=c(points,s)
    }
  }
  if(s<=T) return(points)
  else{
    points=points[points<=T]
    return(points)
  }
}
######## Setting ###########
set.seed(1)
ind<-function(x) ifelse(x<=0.01,1,0)
int0=function(x) 20000*(x+0.001)*exp(1-500*x)*ind(x)
int1=function(x) -15000*(x+0.001)*exp(1-500*x)*ind(x)
onesim<-function(n,tmax,netstr){
  
  points=vector("list",n)
  # level1
  for(m in level1){
    base1<-function(x) alpha[m]*sin(2*5*pi*x/tmax)+alpha[m]
    points[[m]]=gen0(base1,tmax,100)
  } 
  
  # level2
  lambda.m=2000
  for(m in level2){
    points.2=vector()
    s=0
    base<-function(x) alpha[m]*sin(2*5*pi*x/tmax)+alpha[m]
    while(s<tmax){
      u=runif(1)
      w=-log(u)/lambda.m
      s=s+w
      D=runif(1)
      lambda.m0=base(s)
      posedg=which(netstr[,m]==1)
      if(length(posedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(posedg,function(k) sum(int0(s-points[[k]][points[[k]]<s])) ) )
      }
      negedg=which(netstr[,m]==-1)
      if(length(negedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(negedg,function(k) sum(int1(s-points[[k]][points[[k]]<s])) ) )
      }
      
      if(D*lambda.m<=lambda.m0){
        points.2=c(points.2,s)
      }
    }
    points[[m]]=points.2[points.2<=tmax]
  }
  
  # level3
  lambda.m=2000
  
  for(m in level3){
    points.3=vector()
    s=0
    base<-function(x) alpha[m]*sin(2*5*pi*x/tmax)+alpha[m]
    while(s<tmax){
      u=runif(1)
      w=-log(u)/lambda.m
      s=s+w
      D=runif(1)
      lambda.m0=base(s)
      posedg=which(netstr[,m]==1)
      if(length(posedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(posedg,function(k) sum(int0(s-points[[k]][points[[k]]<s])) ) )
      }
      negedg=which(netstr[,m]==-1)
      if(length(negedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(negedg,function(k) sum(int1(s-points[[k]][points[[k]]<s])) ) )
      }
      if(D*lambda.m<=lambda.m0){
        points.3=c(points.3,s)
      }
    }
    points[[m]]=points.3[points.3<=tmax]
  }
  
  # level4
  lambda.m=2000
  
  for(m in level4){
    points.4=vector()
    s=0
    base<-function(x) alpha[m]*sin(2*5*pi*x/tmax)+alpha[m]
    while(s<tmax){
      u=runif(1)
      w=-log(u)/lambda.m
      s=s+w
      D=runif(1)
      lambda.m0=base(s)
      posedg=which(netstr[,m]==1)
      if(length(posedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(posedg,function(k) sum(int0(s-points[[k]][points[[k]]<s])) ) )
      }
      negedg=which(netstr[,m]==-1)
      if(length(negedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(negedg,function(k) sum(int1(s-points[[k]][points[[k]]<s])) ) )
      }
      if(D*lambda.m<=lambda.m0){
        points.4=c(points.4,s)
      }
    }
    points[[m]]=points.4[points.4<=tmax]
  }
  
  # level5
  lambda.m=2000
  
  for(m in level5){
    points.5=vector()
    s=0
    base<-function(x) alpha[m]*sin(2*5*pi*x/tmax)+alpha[m]
    while(s<tmax){
      u=runif(1)
      w=-log(u)/lambda.m
      s=s+w
      D=runif(1)
      lambda.m0=base(s)
      posedg=which(netstr[,m]==1)
      if(length(posedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(posedg,function(k) sum(int0(s-points[[k]][points[[k]]<s])) ) )
      }
      negedg=which(netstr[,m]==-1)
      if(length(negedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(negedg,function(k) sum(int1(s-points[[k]][points[[k]]<s])) ) )
      }
      if(D*lambda.m<=lambda.m0){
        points.5=c(points.5,s)
      }
    }
    points[[m]]=points.5[points.5<=tmax]
  }
  
  # level6
  lambda.m=2000
  
  for(m in level6){
    points.6=vector()
    s=0
    base<-function(x) alpha[m]*sin(2*5*pi*x/tmax)+alpha[m]
    while(s<tmax){
      u=runif(1)
      w=-log(u)/lambda.m
      s=s+w
      D=runif(1)
      lambda.m0=base(s)
      posedg=which(netstr[,m]==1)
      if(length(posedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(posedg,function(k) sum(int0(s-points[[k]][points[[k]]<s])) ) )
      }
      negedg=which(netstr[,m]==-1)
      if(length(negedg)!=0){
        lambda.m0=lambda.m0+sum(sapply(negedg,function(k) sum(int1(s-points[[k]][points[[k]]<s])) ) )
      }
      if(D*lambda.m<=lambda.m0){
        points.6=c(points.6,s)
      }
    }
    points[[m]]=points.6[points.6<=tmax]
  }
  
  points
}

##########################################
########## estimation ####################
##########################################
library(splines)
bs0 <- function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
                 Boundary.knots = range(x)) 
{
  if(length(x)==0) return(matrix(0,1,length(knots)+2+degree-1))
  
  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1) 
    stop("'degree' must be integer >= 1")
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  }
  else FALSE
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - ord + (1L - intercept)
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d", 
                       ord - (1L - intercept)), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                         2L)[-c(1L, nIknots + 2L)]
      quantile(x[!outside], knots)
    }
  }
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if (any(outside)) {
    return(matrix(0,1,length(knots)+2+degree-1))
  }
  else basis <- splineDesign(Aknots, x, ord)
  if (!intercept) 
    basis <- basis[, -1L, drop = FALSE]
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots, 
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
}
int.knots=seq(0,20,length.out=18)
int.knots=int.knots[2:17]
BX0 <- function(x) bs0(x, degree=3, knots=int.knots, intercept = TRUE, Boundary.knots=c(0, 20))

knots=seq(0,0.01,0.00125)
sf<-function(x){
  sf1=ifelse(x<knots[2]& x>=knots[1],1,0)
  sf2=ifelse(x<knots[3]& x>=knots[2],1,0)
  sf3=ifelse(x<knots[4]& x>=knots[3],1,0)
  sf4=ifelse(x<knots[5]& x>=knots[4],1,0)
  sf5=ifelse(x<knots[6]& x>=knots[5],1,0)
  sf6=ifelse(x<knots[7]& x>=knots[6],1,0)
  sf7=ifelse(x<knots[8]& x>=knots[7],1,0)
  sf8=ifelse(x<knots[9]& x>=knots[8],1,0)
  cbind(sf1,sf2,sf3,sf4,sf5,sf6,sf7,sf8)
}
G.est=function(h0,n0,n1){
  x=seq(0,19.9999,0.0001)
  M=length(h0)
  psi0=BX0(x)
  psi=lapply(1:M,function(m){
    t(sapply(1:length(x),function(j) colSums(sf(x[j]-h0[[m]][h0[[m]]<x[j]&h0[[m]]>x[j]-.01]))  ))
  } )
  G1.0=t(psi0)%*%psi0*0.0001
  G1.1=lapply(1:M,function(m){
    t(psi0)%*%psi[[m]]*0.0001
  })
  
  
  G2.1=lapply(1:M,function(m) t(psi[[m]])%*%psi[[m]]*0.0001 )
  
  G2.2=sapply(1:(M-1),function(m1){
    lapply((m1+1):M,function(m2){
      t(psi[[m1]])%*%psi[[m2]]*0.0001
    })
  })
  G=matrix(0,n1*M+n0,n1*M+n0)
  G[1:n0,1:n0]=G1.0
  for(j in 1:M) G[1:n0,(n1*j+n0-n1+1):(n1*j+n0)]=G1.1[[j]]
  
  for(i in 1:M){
    G[(n1*i+n0-n1+1):(n1*i+n0),(n1*i+n0-n1+1):(n1*i+n0)]=G2.1[[i]]
  }
  for(i in 1:(M-1)){
    for(j in (i+1):M){
      G[(n1*i+n0-n1+1):(n1*i+n0),(n1*j+n0-n1+1):(n1*j+n0)]=G2.2[[i]][[j-i]]
    }
  }
  
  
  for(i in 1:(n1*M+n0)){
    for(j in 1:(i-1)){
      G[i,j]=G[j,i]
    }
  }
  G
}

b.est=function(h0,dis,i,n1){
  x=seq(0,dis[1],0.0001)
  M=length(h0)
  b1=colSums(BX0(h0[[i]]))
  b1.0=sapply(1:M,function(m) sapply(1:n1,function(k) sum(sapply(1:length(h0[[i]]),function(j) 
    sum(sf(h0[[i]][j]-h0[[m]][(h0[[m]]<h0[[i]][j])&(h0[[m]]>h0[[i]][[j]]-0.01)])[,k])))))
  
  b=b1
  for(i in 1:M) b=c(b,b1.0[,i])
  b
}

###################################
########## Tunning ################
###################################
grplas0<-function(G,b,lambda,M,n1,n0){
  b.t=list()
  b.t0=b[1:n0]
  for(m in 1:M) b.t[[m]]=b[(n1*m-n1+n0+1):(n1*m+n0)]
  U=lapply(1:M,function(m){
    G0=G[(n1*m-n1+n0+1):(n1*m+n0),(n1*m-n1+n0+1):(n1*m+n0)]
    Q=svd(G0)$u
    D=svd(G0)$d
    Q%*%diag(D^(-1/2))%*%t(Q)
  } )
  G0=G[1:n0,1:n0]
  Q=svd(G0)$u
  D=svd(G0)$d
  u0=Q%*%diag(D^(-1/2))%*%t(Q)
  H0=lapply(1:M,function(m) t(u0)%*%G[1:n0,(n1*m-n1+n0+1):(n1*m+n0)]%*%U[[m]])
  H=lapply(1:M,function(m1){
    lapply(1:M,function(m2) t(U[[m1]])%*%G[(n1*m1-n1+n0+1):(n1*m1+n0),(n1*m2-n1+n0+1):(n1*m2+n0)]%*%U[[m2]] )
  })
  b.t0=u0%*%b.t0
  for(m in 1:M) b.t[[m]]=U[[m]]%*%b.t[[m]]
  
  mu.new=mu=rep(0,n0)
  beta.new=beta=matrix(0,M,n1)
  #z.new=z
  delta=10
  while(delta>10^(-6)){
    other0=lapply(1:M,function(m) H0[[m]]%*%beta[m,]  )
    other0.0=other0[[1]]
    for(m in 2:M) other0.0=other0.0+other0[[m]]
    mu.new=b.t0-other0.0
    for(m in 1:M){
      other=lapply(1:M,function(m1) H[[m]][[m1]]%*%beta.new[m1,])
      other.0=-other[[m]]
      for(m1 in 1:M) other.0=other.0+other[[m1]]
      r=b.t[[m]]-t(H0[[m]])%*%mu.new-other.0
      r.norm=sqrt(sum(r^2))
      beta.new[m,]=if(r.norm>lambda) (r.norm-lambda)/r.norm*r else 0
    }
    delta=(sum((mu.new-mu)^2)+sum((beta.new-beta)^2))/(sum(mu.new^2)+sum(beta.new^2))
    mu=mu.new
    beta=beta.new
  }
  
  a.0=u0%*%mu
  a.1=sapply(1:M,function(m) U[[m]]%*%beta[m,])
  a1=a.0
  for(m in 1:M) a1=c(a1,a.1[,m])
  a1
}
grplas<-function(G,b,lambda,M,n1){
  z=list()
  z0=b[1]
  for(m in 1:M) z[[m]]=b[(n1*m-(n1-2)):(n1*m+1)]
  U=lapply(1:M,function(m){
    G0=G[(n1*m-(n1-2)):(n1*m+1),(n1*m-(n1-2)):(n1*m+1)]
    Q=svd(G0)$u
    D=svd(G0)$d
    Q%*%diag(D^(-1/2))%*%t(Q)
  } )
  u0=G[1,1]^(-1/2)
  z0=z0*u0
  for(m in 1:M) z[[m]]=t(z[[m]]%*%U[[m]])
  
  mu.new=mu=0
  beta.new=beta=matrix(0,M,n1)
  z.new=z
  delta=10
  while(delta>10^(-6)){
    z0.new=z0-sum(sapply(1:M,function(m) u0*G[1,(n1*m-(n1-2)):(n1*m+1)]%*%U[[m]]%*%beta[m,]  ))
    mu.new=z0.new
    for(m in 1:M){
      other=lapply(1:M,function(m1) t(U[[m]])%*%G[(n1*m-(n1-2)):(n1*m+1),(n1*m1-(n1-2)):(n1*m1+1)]%*%U[[m1]]%*%beta.new[m1,])
      other.0=-other[[m]]
      for(m1 in 1:M) other.0=other.0+other[[m1]]
      z.new[[m]]=z[[m]]-t(U[[m]])%*%G[(n1*m-(n1-2)):(n1*m+1),1]%*%mu.new*u0-other.0
      z.norm=sqrt(sum(z.new[[m]]^2))
      beta.new[m,]=if(z.norm>lambda) (z.norm-lambda)/z.norm*z.new[[m]] else 0
    }
    delta=((mu.new-mu)^2+sum((beta.new-beta)^2))/(mu.new^2+sum(beta.new^2))
    mu=mu.new
    beta=beta.new
  }
  
  a.0=mu*u0
  a.1=sapply(1:M,function(m) U[[m]]%*%beta[m,])
  a1=a.0
  for(m in 1:M) a1=c(a1,a.1[,m])
  a1
}
econ3<-function(x,netstr){
  edge.id=which(netstr!=0)
  npos=length(edge.id)
  nneg=length(netstr)-npos
  
  tedges=sum(x!=0)
  tpos=sum(x[edge.id]==netstr[edge.id])
  tneg=tedges-tpos
  fnr=(npos-tpos)/npos
  fpr=(tneg)/nneg
  f1=2/(npos/tpos+tedges/tpos)
  tpc=matrix(0,2,2)
  tpc[1,]=c(tpos,tneg)
  tpc[2,]=c(npos-tpos,nneg-tneg)
  list(total.egdes=tedges,tpc=tpc,tfnr=fnr,fpr=fpr,f1=f1)
}