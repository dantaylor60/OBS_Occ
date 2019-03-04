require(truncnorm)
require(Rcpp)
require(RcppArmadillo)
require(combinat)
sourceCpp("CppFns.cpp")

#Model space enumeration functions
.allterm.generator <-function(v,r){
  do.call(rbind,lapply(1:r,FUN=function(v,x){
    t(xsimplex(p=v,n=x))},v=v))
}
.nodes <- function(x,O1.names){
  if(sum(x)==0){
    c(nam="Int",form.alias="(Intercept)")
  }else{
    nam = O1.names[x!=0]
    nam.alias = character(0)
    xin = x[x!=0]
    for(i in 1:length(nam)){
      nam[i] = ifelse(xin[i]==1,nam[i],paste(nam[i],"^",xin[i],sep=""))  
      nam.alias[i] = ifelse(xin[i]==1,nam[i],paste("I(",nam[i],")",sep=""))  
    }
    if(length(nam>1)){
      nam=paste(nam,collapse="*")
      nam.alias=paste(nam.alias,collapse=":")
      c(nam=nam,form.alias=nam.alias)
    }else{
      c(nam=nam,form.alias=nam.alias)
    }
  }
}
.node.edges <- function(x,polymat,vertices,degree.vec,degree.max){
  if(degree.vec[x]< degree.max){
    do.call(rbind,lapply(1:dim(polymat)[2],FUN=function(u){
      pos.x = which(vertices%in%x)
      poly.xto = polymat[pos.x,]
      poly.xto[u] = poly.xto[u]+1
      c(x,names(which(apply(polymat,MARGIN=1,identical,y=poly.xto))))
    }))
  }
}
.MF.DAG <- function(X.linear,maxdeg,exclusions=NA,MF.formula=NA){
  deg1names = names(X.linear)
  var.class = c(lapply(1:length(deg1names),function(x)class(X.linear[,x])),recursive=T)
  categorical = deg1names[var.class=="factor"]
  if(!identical(categorical,character(0)) && maxdeg>1){ 
    categorical = paste(categorical,"^2",sep="")
    exclusions = c(categorical,exclusions)
    exclusions = unique(exclusions[!is.na(exclusions)])
  }
  nvars = length(deg1names)
  polyterms = rbind(rep(0,nvars),.allterm.generator(v=nvars,r=maxdeg))
  vertices = apply(polyterms,MARGIN=1,FUN = .nodes, O1.names=deg1names )
  form.names = vertices["form.alias",]
  vertices = vertices["nam",]
  rownames(polyterms) = vertices 
  degree = rowSums(polyterms)
  edges = do.call(rbind,lapply(vertices,FUN = .node.edges,polymat=polyterms,
                               vertices=vertices,degree.vec=degree,degree.max=maxdeg))
  if(sum(is.na(exclusions))==0){
    
    mindeg.excl = min(degree[vertices%in%exclusions])
    
    for(i in mindeg.excl:(maxdeg+1)){
      rmv = do.call(rbind,lapply(exclusions,FUN=.node.edges,
                                 polyterms,vertices,degree,maxdeg))
      exclusions = unique(c(exclusions,rmv[,2]))
    }
    form.names = form.names[!(vertices%in%exclusions)]
    degree = degree[!(vertices%in%exclusions)]
    vertices = vertices[!(vertices%in%exclusions)]
    edges = edges[!(edges[,2]%in%exclusions),]
  }
  if(class(MF.formula)=="formula"){
    nn1 = gsub(pattern=" ",replacement="",x=as.character(MF.formula)[3],fixed=T)
    vert.MF = strsplit(nn1,split="+",fixed=T)[[1]]
    vert.MF = unique(c("Int",vert.MF))
    form.names = form.names[(vertices%in%vert.MF)]
    edges = edges[(edges[,2]%in%vert.MF),]
    degree = degree[(vertices%in%vert.MF)]
    vertices = vertices[(vertices%in%vert.MF)]
  }
  
  parents = lapply(vertices,FUN=function(x){
    c(edges[edges[,2]==x,1])
  })  
  names(parents)=vertices
  
  return(list(edges=edges,vertices=vertices,degree=degree,
              parents=parents,form.names=form.names))
}
.model.enumerate<-function(MB,MF){
  
  bit.increase<-function(x){
    if(prod(x)==TRUE){return(c(rep(FALSE,length(x)),TRUE))}
    index<-1
    flag<-1
    while(flag){
      if(x[index]==FALSE){x[index]<-TRUE;flag<-0}else{x[index]<-FALSE;index<-index+1}
    }
    return(x)
  }
  
  v<-MF$vertices
  v.base<-MB$vertices
  v.test<-MF$vertices[!(v%in%v.base)]
  degree<-MF$degree
  model.expand<-function(model,D){
    vert.in<-c(v.base,v.test[model])
    points.to<-unique(MF$edges[MF$edges[,1]%in%vert.in,2])
    points.to<-points.to[!(points.to%in%vert.in)]
    #print(points.to)
    #print(vert.in)
    v.loc<-NULL
    for(pt in points.to){
      if((prod(MF$parents[[pt]]%in%vert.in)==1)){
        v.loc<-c(v.loc,pt)
      }
    } 
    #print(v.loc)
    L<-length(v.loc)
    degree.loc<-rep(NA,L)
    for(i in 1:L){
      degree.loc[i]<-degree[v==v.loc[i]]
    }
    #print(vert.in)
    #print(v.loc)
    #print(degree.loc)
    v.loc<-v.loc[degree.loc>D]
    degree.loc<-degree.loc[degree.loc>D]
    L<-length(v.loc)
    if(L==0){return(model)}
    models<-model
    for(d in min(degree.loc):max(degree.loc)){
      #print(d)
      v.loc2<-v.loc[degree.loc==d]
      L<-length(v.loc2)
      x<-rep(0,L)
      for(i in 1:(2^L-1)){
        x<-bit.increase(x)
        #print(x)
        model.new<-model
        model.new[v.test%in%v.loc2]<-as.logical(x)
        models.new<-model.expand(model.new,D=d)
        models<-rbind(models,models.new)
      }
    }
    return(models)
  }
  model<-rep(FALSE,length(v.test))
  models<-model.expand(model,D=0)
  rownames(models)<-NULL
  colnames(models)<-v.test
  return(models)
}

#occupancy model fitting functions
.formatData <- function(data,form,Ji){
  if(form!=as.formula(. ~ 1)){
    facts = names(data)[sapply(data,class)=="factor"]
    data[!(names(data)%in%facts)] = scale(data[!(names(data)%in%facts)])
    vvz = character(0)
    namvars = strsplit(as.character(form)[3],split=" + ",fixed=T)[[1]]
    for(j in namvars){
      comps=strsplit(j,split=" * ",fixed=T)[[1]]
      if(length(comps)==1){
        if(comps%in%facts){
          vvz = c(vvz,j)
        }else{
          vvz = c(vvz,paste0("I(",j,")"))
        }
      }else{
        if(sum(comps%in%facts)==0){
          vvz = c(vvz,paste0("I(",comps,")",collapse=":"))
        }else{
          vvz = c(vvz,paste(unlist(lapply(comps,function(xx){
            if(xx%in%facts){
              xx
            }else{
              paste0("I(",xx,")")
            }
          })),collapse=":"))
        }
      }
    }
    ff = as.formula(paste0("~",paste(vvz,collapse="+")))
  }else{
    ff = as.formula(~1)
  }
  return(model.matrix(terms(ff,keep.order=T),data))
}
.negLLwrap = function(Y,X,Q){
  Ji <- apply(Y,1,function(x){sum(!is.na(x))})
  Jtot <- sum(Ji)
  yvec <-c(na.exclude(c(t(Y))))
  
  function(param) {
    n=nrow(X)
    Ji
    betaz = matrix(param[1:ncol(X)],ncol=1)
    psi = pnorm(X %*% betaz)
    
    betay = matrix(param[-(1:ncol(X))],ncol=1)
    p = pnorm(Q%*%betay)
    
    logL = rep(NA,n)
    for (i in 1:n) {
      s=(sum(Ji[0:(i-1)],1):sum(Ji[1:i]))
      terms = p[s]^yvec[s] * (1-p[s])^(1-yvec[s])
      logL[i] = log( psi[i] * prod(terms)  + ifelse(sum(yvec[s])==0,1,0)*(1-psi[i])  )
    }
    return((-1)*sum(logL))
  }
}
.OccuFit.wrap <- function(dataZ,dataY,Y,pxda=T,nams.MBz,nams.MBy,
                          Storage,PostInfo){
  #-------------define constant objects
  N <- dim(dataZ)[1]
  Ji <- apply(Y,1,function(x){sum(!is.na(x))})
  Jtot <- sum(Ji)
  yvec <-c(na.exclude(c(t(Y))))
  z.ini = as.numeric(rowSums(Y,na.rm=T)>0)
  ind.zu = which(z.ini==0)
  
  cnst.z = ncol(dataZ)/(2*N)
  X0 = dataZ[,nams.MBz,drop=F]
  Pz0 = X0%*%solve(crossprod(X0))%*%t(X0)
  SigmaInv.z = t(dataZ) %*% ( diag((1+cnst.z),N) - cnst.z*Pz0) %*% dataZ
  Sigma.z = solve(SigmaInv.z)
  Hz = tcrossprod(Sigma.z,dataZ)
  PostInfo$Sigma.z = Sigma.z
  
  cnst.y = ncol(dataY)/(2*Jtot)
  Q0 = dataY[,nams.MBy,drop=F]
  Py0 = Q0%*%solve(crossprod(Q0))%*%t(Q0)
  SigmaInv.y = t(dataY) %*% ( diag((1+cnst.y),Jtot) - cnst.y*Py0) %*% dataY
  Sigma.y = solve(SigmaInv.y)
  Hy = tcrossprod(Sigma.y,dataY)
  PostInfo$Sigma.y = Sigma.y
  
  #-------------define intermediate functions
  #--- pxda function
  gsqsamplewrap <- function(pxda){
    function(resp,mu,n){
      if(pxda){
        sqrt(rgamma(1,shape=n/2,rate=sum((resp-mu)^2)/2))
      }else{
        1
      }
    }
  }
  gsample<-gsqsamplewrap(pxda=pxda)
  #--- psi.star function
  wrap.psi <- function(Ji,vecposz){
    function(py,psi){
      unlist(lapply(vecposz,function(posz){
        s=sum(Ji[0:(posz-1)],1)
        return(psi[posz]*prod(1-py[s:(s+Ji[posz]-1)])/( psi[posz]*prod(1-py[s:(s+Ji[posz]-1)]) + (1-psi[posz])))
      }))}
  }
  fn.psistar <- wrap.psi(Ji,vecposz=ind.zu)
  
  function(iter){
    z.t <- Storage$zsamples[iter,]
    betaz.t <- matrix(Storage$betaz[iter,],ncol=1)
    betay.t <- matrix(Storage$betay[iter,],ncol=1)
    
    muY <- as.vector(dataY%*%betay.t)
    w.t <- rtruncnorm(Jtot, a=log(yvec*rep(z.t,Ji)), b=-log((1-yvec)*rep(z.t,Ji)), 
                      mean = muY, sd = 1)
    
    gW <- gsample(resp=w.t,mu=muY,n=Jtot)
    betayRB.t <- Hy%*%(gW*w.t)
    
    
    muZ <- as.vector(dataZ%*%betaz.t)
    p.t <- pnorm(muY)
    psi.t <- pnorm(muZ)
    psi.star <- fn.psistar(py=p.t,psi=psi.t)
    z.t[ind.zu] = rbinom(length(ind.zu),1,psi.star)  
    v.t <- rtruncnorm(N, a=log(z.t), b=-log(1-z.t), mean = muZ, sd = 1)
    gV <- gsample(resp=v.t,mu=muZ,n=N)
    betazRB.t <- Hz%*%(gV*v.t)
    
    #-----------------------------------------------
    Storage$zsamples[(iter+1),ind.zu] = z.t[ind.zu]
    Storage$betaz[(iter+1),] = rmvnormArma(1, betazRB.t, Sigma.z)
    Storage$betay[(iter+1),] = rmvnormArma(1, betayRB.t, Sigma.y)
    
    PostInfo$wsamples[iter,] <- w.t
    PostInfo$vsamples[iter,] <- v.t
  }
}

OccuFit.fn <- function(Ytr,dataZ.tr,dataY.tr,formZ,formY,niter,pxda,MBz="(Intercept)",MBy="(Intercept)"){
  
  Xtr<-.formatData(data=dataZ.tr,form=formZ)
  Qtr<-.formatData(data=dataY.tr,form=formY)
  
  #-------------define storage objects
  StoreEstim <- setRefClass("storest",fields=list(zsamples="matrix",betaz="matrix",betay="matrix"))
  StorePosterior <- setRefClass("storelat",fields=list(vsamples="matrix",wsamples="matrix",
                                                       Sigma.z="matrix",Sigma.y="matrix"))
  
  StorageEst <- StoreEstim$new()
  StorageEst$zsamples = matrix(1,ncol=nrow(Ytr),nrow=(niter+1))
  StorageEst$betaz = matrix(0,ncol=ncol(Xtr),nrow=(niter+1))
  StorageEst$betay = matrix(0,ncol=ncol(Qtr),nrow=(niter+1))
  
  StoragePost <- StorePosterior$new()
  StoragePost$vsamples = matrix(0,ncol=nrow(Xtr),nrow=(niter))
  StoragePost$wsamples = matrix(0,ncol=nrow(Qtr),nrow=(niter))
  StoragePost$Sigma.z = matrix(0,nrow=ncol(Xtr),ncol=ncol(Xtr))
  StoragePost$Sigma.y = matrix(0,nrow=ncol(Qtr),ncol=ncol(Qtr))
  
  #-------------initialize parameters
  negloglikfn <- .negLLwrap(Y=Ytr,X=Xtr,Q=Qtr)
  Xtr=as.matrix(Xtr)
  betazGuess = rep(0,ncol(Xtr))
  betayGuess = rep(0,ncol(Qtr))
  paramGuess = c(betazGuess, betayGuess)
  fitMLE = optim(par=paramGuess, fn=negloglikfn, 'SANN',control=list(maxit=1000))$par
  StorageEst$betaz[1,] = fitMLE[1:ncol(Xtr)]
  StorageEst$betay[1,] = fitMLE[-(1:ncol(Xtr))]
  z.ini = as.numeric(rowSums(Ytr,na.rm=T)>0)
  
  if(ncol(Xtr)>1){
    StorageEst$zsamples[1,(z.ini==0)] = rbinom(sum(z.ini==0),1,pnorm(Xtr[(z.ini==0),]%*%
                                                                       fitMLE[1:ncol(Xtr)]))
  }else{
    StorageEst$zsamples[1,(z.ini==0)] = rbinom(sum(z.ini==0),1,pnorm(Xtr[(z.ini==0),]*fitMLE))
  }

  #---------------initialize estimation function
  fitoccufn <- .OccuFit.wrap(dataZ=Xtr,dataY=Qtr,Y=Ytr,pxda=pxda,nams.MBz=MBz,
                             nams.MBy=MBy,Storage=StorageEst,PostInfo=StoragePost)
  
  #---------------Gibbs sampler loop
  pb <- txtProgressBar(min = 0, max = niter, style = 3)
  for(i in 1:niter){
    fitoccufn(iter=i)
    setTxtProgressBar(pb, i)
  }
  return(list(out.est=StorageEst,out.post=StoragePost,X=Xtr,Q=Qtr))
}

getmarginal.Y <- function(mcmc.outpost,Y,X,Q,betaz,betay,burnin,thin){
  require(mvtnorm)
  
  n <- nrow(X)
  nsim <- nrow(mcmc.outpost$vsamples)
  Jtot <- nrow(Q)
  pz <- ncol(X)
  py <- ncol(Q)
  
  insamples <- seq(burnin,nsim,by=thin)
  #calculate parameter posterior means
  mupost.betaz <- apply(mcmc.outpost$vsamples[insamples,],1,function(x){
    tcrossprod(mcmc.outpost$Sigma.z,X)%*%matrix(x,ncol=1)
  })
  mupost.betay <- apply(mcmc.outpost$wsamples[insamples,],1,function(x){
    tcrossprod(mcmc.outpost$Sigma.y,Q)%*%matrix(x,ncol=1)
  })
  
  #calculate parameter priors and posteriors for specific parameter values
  if(pz==1){
    pr.parsz <- 0
    post.parsz <- unlist(lapply(mupost.betaz,function(mu)dnorm(betaz,mean=mu,sd=sqrt(mcmc.outpost$Sigma.z[1]),log=T)))
  }else{
    cnst.z = pz/(2*n)
    X0t = X[,1,drop=F]
    Sigma.prz <- solve(0.5*cnst.z*t(X)%*%(diag(n)-(1/n)*cnst.z*tcrossprod(X0t))%*%X)
    pr.parsz <- dmvnorm(betaz, sigma = Sigma.prz, log = TRUE)
    post.parsz <- apply(mupost.betaz,2,function(mu)dmvnorm(x=betaz,mean=mu,sigma=mcmc.outpost$Sigma.z,log=T))
  }
  
  if(py==1){
    pr.parsy <- 0
    post.parsy <- unlist(lapply(mupost.betay,function(mu)dnorm(betay,mean=mu,sd=sqrt(mcmc.outpost$Sigma.y[1]),log=T)))
  }else{
    cnst.y = py/(2*n)
    Q0t = Q[,1,drop=F]
    Sigma.pry <- solve(0.5*cnst.y*t(Q)%*%(diag(Jtot)-(1/Jtot)*cnst.y*tcrossprod(Q0t))%*%Q)
    pr.parsy <- dmvnorm(betay, sigma = Sigma.pry, log = TRUE)
    post.parsy <- apply(mupost.betay,2,function(mu)dmvnorm(x=betay,mean=mu,sigma=mcmc.outpost$Sigma.y,log=T))
  }
  
  nll <- .negLLwrap(Y=Y,X=X,Q=Q)
  llik <- (-nll(c(betaz,betay)))
  logprior <- (pr.parsz+pr.parsy)
  maxlog.zy <- max(post.parsz) + max(post.parsy)
  logpost <- log(mean(exp(post.parsz+post.parsy-maxlog.zy)))+maxlog.zy
  
  return(llik+logprior-logpost)
  
}


#load Mallard Y, site and site-obs data, as well as model spaces for Z and Y
load("MallardInfo.RData")
for(i in 1:length(list.Mallard)){
  assign(names(list.Mallard)[i],list.Mallard[[i]])
}
rm("list.Mallard","i")

#start the 
nsim=10000
test=OccuFit.fn(Ytr=Y,dataZ.tr=X.linear,dataY.tr=Q.linear,formZ=as.formula(vec.Mspz[50]),
                formY=as.formula(vec.Mspy[4]),niter=nsim,pxda=F,MBz="(Intercept)",MBy="(Intercept)")
logmargY <- getmarginal.Y(mcmc.outpost = test$out.post, Y = Y, X = test$X, Q = test$Q,
                          betaz = test$out.est$betaz[100,], betay = test$out.est$betay[100,],
                          burnin=round(0.2*nsim),thin=5)