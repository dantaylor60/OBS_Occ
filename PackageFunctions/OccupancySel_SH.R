require(glmnet)
require(mvtnorm)
require(truncnorm)
require(Rcpp)
require(RcppArmadillo)
require(RcppEigen)
require(inline)
require(fastmatch)
sourceCpp("CppFns.cpp")
#library(package="UtilModSel1")
source("wrapperStrongHeredity.R")


Occup.IPSel <- function(Y,full.formulaZ,dataZ,base.formulaZ=as.formula(.~ 1),full.formulaY,dataY,base.formulaY=as.formula(.~ 1),
                        SH=TRUE,maxdegZ=0,maxdegY=0,nodes.to.removeZ=NULL,nodes.to.removeY=NULL,model.prior.type="EPP",pxda=F,
                        model.prior.pars="children",beta.prior.type="IP",beta.prior.pars=list(alpha=1,nu=1),eff=TRUE,centerscale=TRUE,
                        start.formulaZ=as.formula(.~ 1),start.formulaY=as.formula(.~ 1),niter=5000,ordy.by.site=T,
                        search.start=T,pshuffle=0.1,burnin=1000,thin=5){
  #-------------------------------------
  # dataZ and dataY must include response in 1st column and predictors in the rest
  # dataY must be ordered ordered first by site then by survey if ordy.by.site=T
  # (i.e. c(site1.surv1, site1.surv2, ..., site1.survJ1, site2.surv1,site2.surv2,..., siten.survJn))
  #
  #
  #---------------------------------------------
  #Define mutable objects for Y and for Z (ModelStorage, ModelSequence, logPostVectorm, MatParRB, LatentStore)
  #---------------------------------------------
  ModStore <- setRefClass("ModStore",fields=list(model="list"))
  ModelStorageZ <- ModStore$new()
  ModelStorageY <- ModStore$new()
  #---------------------------------------------
  ModSeq <- setRefClass("ModSeq",fields=list(modseq="vector"))
  model.sequenceZ <- ModSeq$new() 
  model.sequenceY <- ModSeq$new()
  model.sequenceZ$modseq <- vector(length=niter)
  model.sequenceY$modseq <- vector(length=niter)
  #---------------------------------------------
  logPostVector <- setRefClass("logPost",fields=list(logpost="list"))
  logPostZ <- logPostVector$new()
  logPostY <- logPostVector$new()
  #---------------------------------------------
  MatParRB <- setRefClass("MatParRB",fields=list(ParRB="matrix"))
  MatParRBZ <- MatParRB$new()
  MatParRBY <- MatParRB$new()
  #---------------------------------------------
  #Initialize objects
  #---------------------------------------------
  N <- dim(dataZ)[1]
  Ji <- apply(Y,1,function(x){sum(!is.na(x))})
  Jtot <- sum(Ji)
  cellid <- 1:N
  cellvec <- rep(cellid,Ji)
  yvec <-c(na.exclude(c(t(Y))))
  
  if(!ordy.by.site){
    ordmat=matrix(NA,ncol=ncol(Y),nrow=nrow(Y))
    ordmat[!is.na(Y)]=1:length(yvec)
    ord = na.exclude(c(t(ordmat)))
    dataY=dataY[ord,]
    rm(list="ordmat")
  }
  
  z.t = as.numeric(rowSums(Y,na.rm=T)>0)
  ind.zu = which(z.t==0)
  z.t[ind.zu] = rbinom(length(ind.zu),1,0.5)
    
  #---------------Additional mutable objects for latents
  LatStore <- setRefClass("LatStore",fields=list(latents="matrix"))
  StoreW <- LatStore$new()
  StoreV <- LatStore$new()
  seqStorelat <- seq((burnin+1),niter,by=thin)
  lseq <- length(seqStorelat)
  StoreV$latents <- matrix(NA,ncol=N,nrow=length(seqStorelat))
  StoreW$latents <- matrix(NA,ncol=Jtot,nrow=length(seqStorelat))
  ind.store <- 1
  
  #------  PXDA function
  
  gsqsamplewrap <- function(pxda){
    if(pxda){
      function(resp,mu,n) sqrt(rgamma(1,shape=n/2,rate=sum((resp-mu)^2)/2))
      #function(a,b,n){
      #  #a=y'y; b = y'X beta.hat; beta.hat=(X'X+Omega)^{-1} X'y
      #    sqrt(rgamma(1,shape=n/2,rate=(a-b)/2))
    }else{
        1
    }
  }
  
  gsample<-gsqsamplewrap(pxda=pxda)
  
  v.t <- rtruncnorm(N, a=log(z.t), b=-log(1-z.t), mean = 0, sd = 1)
  w.t <- rtruncnorm(Jtot, a=log(yvec*rep(z.t,Ji)), b=-log((1-yvec)*rep(z.t,Ji)), 
                    mean = 0, sd = 1)
  functionselZ <- psr(full.formula=full.formulaZ,data=dataZ,base.formula=base.formulaZ,SH=SH,maxdeg=maxdegZ,
                    nodes.to.remove=nodes.to.removeZ,model.prior.type=model.prior.type,
                    model.prior.pars=model.prior.pars,beta.prior.type=beta.prior.type,beta.prior.pars=beta.prior.pars,
                    eff=eff,centerscale=centerscale,start.formula=start.formulaZ,niter=niter,ModelStorage=ModelStorageZ,
                    model.sequence=model.sequenceZ,logpostvec=logPostZ,MatParRB=MatParRBZ,latrespini=v.t,shuffle=pshuffle)
  
  functionselY <- psr(full.formula=full.formulaY,data=dataY,base.formula=base.formulaY,SH=SH,maxdeg=maxdegY,
                    nodes.to.remove=nodes.to.removeY,model.prior.type=model.prior.type,
                    model.prior.pars=model.prior.pars,beta.prior.type=beta.prior.type,beta.prior.pars=beta.prior.pars,
                    eff=eff,centerscale=centerscale,start.formula=start.formulaY,niter=niter,ModelStorage=ModelStorageY,
                    model.sequence=model.sequenceY,logpostvec=logPostY,MatParRB=MatParRBY,latrespini=w.t,shuffle=pshuffle)
  
  Model.currentY <- as.list(environment(functionselY))$Model.current
  nodesY <- as.list(environment(functionselY))$nodes
  matFullY <- cbind(as.list(environment(functionselY))$matM0,as.list(environment(functionselY))$matMF)
  muY <- as.vector(matFullY%*%Model.currentY$pars)
  p.t = pnorm(muY)
  
  Model.currentZ <- as.list(environment(functionselZ))$Model.current
  nodesZ <- as.list(environment(functionselZ))$nodes
  matFullZ <- cbind(as.list(environment(functionselZ))$matM0,as.list(environment(functionselZ))$matMF)
  muZ <- as.vector(matFullZ%*%Model.currentZ$pars)
  psi.t <- pnorm(muZ)
  psi.star <- rep(1,N)
  wrap.psi <- function(Ji,vecposz){
    function(py,psi){
      unlist(lapply(vecposz,function(posz){
        s=sum(Ji[0:(posz-1)],1)
        return(psi[posz]*prod(1-py[s:(s+Ji[posz]-1)])/( psi[posz]*prod(1-py[s:(s+Ji[posz]-1)]) + (1-psi[posz])))
    }))}
  }
  calc.psi.star <- wrap.psi(Ji,vecposz=ind.zu)
  psi.star[ind.zu] <- calc.psi.star(py=p.t,psi=psi.t)
  z.t[ind.zu] = rbinom(length(ind.zu),1,psi.star[ind.zu])  
  zt.cumvec <- vector(length=N)
  gW <- gV <- 1.0
  jointpostnum <- vector(length=niter)
  pb <- txtProgressBar(min = 0, max = niter, style = 3)
  activeY <- activeZ <-  bY <-  bZ <-  numeric(0)
  for(k in 1:niter){
    #activeY <- which(Model.currentY$pars!=0)
    #bY <- as.numeric(crossprod(w.t,matFullY[,activeY,drop=F])%*%Model.currentY$parsRB[activeY])
    #print(c(sw=sum(w.t^2),bY=bY))
    #gW <- gsample(a=sum(w.t^2),b=Model.currentY$bcst,n=sum(Ji))#b=bY,n=sum(Ji))#
    gW <- gsample(resp=w.t,mu=muY,n=sum(Ji))
    
    Model.currentY <- functionselY(k=k,Model.currentY,lat.response=(gW*w.t))
    muY <- as.vector(matFullY%*%Model.currentY$pars)
    p.t = pnorm(muY)
    w.t <- rtruncnorm(Jtot, a=log(yvec*rep(z.t,Ji)), b=-log((1-yvec)*rep(z.t,Ji)), 
                      mean = muY, sd = 1)
    
    
    #activeZ <- which(Model.currentZ$pars!=0)
    #bZ <- as.numeric(crossprod(v.t,matFullZ[,activeZ,drop=F])%*%Model.currentZ$parsRB[activeZ])
    #print(c(sv=sum(v.t^2),bZ=bZ))
    #gV <- gsample(a=sum(v.t^2),b=Model.currentZ$bcst,n=N)#,b=bZ,n=N)#
    gV <- gsample(resp=v.t,mu=muZ,n=N)
    Model.currentZ <- functionselZ(k=k,Model.currentZ,lat.response=(gV*v.t))
    muZ <- as.vector(matFullZ%*%Model.currentZ$pars)
    psi.t <- pnorm(muZ)
    psi.star[ind.zu] <- calc.psi.star(py=p.t,psi=psi.t)
    z.t[ind.zu] = rbinom(length(ind.zu),1,psi.star[ind.zu])  
    zt.cumvec =  zt.cumvec+z.t
    v.t <- rtruncnorm(N, a=log(z.t), b=-log(1-z.t), mean = muZ, sd = 1)
    
    ind.store <- ind.store+as.numeric((k>seqStorelat[ind.store])*(k<=seqStorelat[lseq]))
    StoreW$latents[ind.store,] <- w.t
    StoreV$latents[ind.store,] <- v.t
    
    jointpostnum[k] <- sum(Model.currentZ$loglik,Model.currentZ$modpr[1],
                           Model.currentY$loglik,Model.currentY$modpr[1])
    setTxtProgressBar(pb, k)
  }
  
  return(list(MF.matZ=as.list(environment(functionselZ))$MF.mat,
              MF.matY=as.list(environment(functionselY))$MF.mat,
              nodesZ=nodesZ,
              nodesY=nodesY,
              zt.cumvec=zt.cumvec/niter,
              model.sequenceZ=model.sequenceZ,
              model.sequenceY=model.sequenceY,
              logPostZ=logPostZ,
              logPostY=logPostY,
              logpostJoint=jointpostnum,
              MatParRBY=MatParRBY,
              MatParRBZ=MatParRBZ,
              ModelStorageY=ModelStorageY,
              ModelStorageZ=ModelStorageZ,
              StoreV=StoreV,
              StoreW=StoreW
              ))
}

