PostIncProb <- function(listmodelnodes,probs,nodes){
  
  MPIPs <- colSums(do.call(rbind,lapply(1:length(probs),function(mod){
    (1:length(nodes)%in%listmodelnodes[[mod]])*probs[mod]
    })))
  names(MPIPs) <- unlist(lapply(nodes,"[[","name"))
  return(MPIPs)
}

summary.fun <- function(objVarSel,priortype=NULL,erg.fn=median){
  if(!is.null(priortype)){
    pr.modZ <- with(objVarSel,unlist(lapply(ModelStorageZ$model,function(mod)mod$modpr[priortype])) )
    pr.modZ0 <- with(objVarSel,unlist(lapply(ModelStorageZ$model,function(mod)mod$modpr[1])) )
    logpostZ <- lapply(1:length(objVarSel$logPostZ$logpost),function(mod){
      objVarSel$logPostZ$logpost[[mod]]-pr.modZ0[mod]+pr.modZ[mod]
    } )
    
    pr.modY <- with(objVarSel,unlist(lapply(ModelStorageY$model,function(mod)mod$modpr[priortype])) )
    pr.modY0 <- with(objVarSel,unlist(lapply(ModelStorageY$model,function(mod)mod$modpr[1])) )
    logpostY <- lapply(1:length(objVarSel$logPostY$logpost),function(mod){
      objVarSel$logPostY$logpost[[mod]]-pr.modY0[mod]+pr.modY[mod]
    } )
  }else{
    logpostZ <- objVarSel$logPostZ$logpost
    logpostY <- objVarSel$logPostY$logpost
  }
  
  models.z <- unlist(with(objVarSel,lapply(lapply(ModelStorageZ$model,"[[","nodes"),
                                           function(mod)paste(unlist(lapply(nodesZ[mod],"[[","name")),collapse=","))))
  models.y <- unlist(with(objVarSel,lapply(lapply(ModelStorageY$model,"[[","nodes"),
                                           function(mod)paste(unlist(lapply(nodesY[mod],"[[","name")),collapse=","))))
  
  names(models.z)=names(models.y)=NULL
  postZp <- unlist(lapply(logpostZ,erg.fn))
  postZp = exp(postZp-max(postZp))
  postZp = postZp/sum(postZp)
  postYp <- unlist(lapply(logpostY,erg.fn))
  postYp = exp(postYp-max(postYp))
  postYp = postYp/sum(postYp)
  postZf <- table(objVarSel$model.sequenceZ$modseq)/length(objVarSel$model.sequenceZ$modseq)
  postYf <- table(objVarSel$model.sequenceY$modseq)/length(objVarSel$model.sequenceY$modseq)

  orderzf=order(postZf,decreasing=T)  
  orderyf=order(postYf,decreasing=T)
  
  orderz=order(postZp,decreasing=T)
  postZp <- data.frame(models.z,post=postZp)[orderz,]
  
  ordery=order(postYp,decreasing=T)
  postYp <- data.frame(models.y,post=postYp)[ordery,]
  
  niter <- length(objVarSel$model.sequenceZ$modseq)
  
  JointMod<- with(objVarSel,cbind(model.sequenceZ$modseq,model.sequenceY$modseq))
  PostFJoint <- table(unlist(apply(JointMod,1,function(x){
    namsZ = as.numeric(strsplit(x[1],split=",")[[1]])
    namsY = as.numeric(strsplit(x[2],split=",")[[1]])
    with(objVarSel,
         paste(c(paste0(unlist(lapply(nodesZ[namsZ],"[[","name")),collapse=",")),
               paste0(unlist(lapply(nodesY[namsY],"[[","name")),collapse=","),sep="_") )
  })))/niter
  PostFJoint <- PostFJoint[order(PostFJoint,decreasing=T)]
  
  #----- post joint P
  JointModP <- with(objVarSel,data.frame(ModZ=model.sequenceZ$modseq[-1],ModY=model.sequenceY$modseq[-1],
                                         logpost=logpostJoint))
  
  JointModP$ModNam <- unlist(apply(JointModP,1,function(x){
    namsZ = as.numeric(strsplit(x[1],split=",")[[1]])
    namsY = as.numeric(strsplit(x[2],split=",")[[1]])
    with(objVarSel,
         paste(paste0("(Z: ",paste0(unlist(lapply(nodesZ[namsZ],"[[","name")),collapse=","),")"),
               paste0("(Y: ",paste0(unlist(lapply(nodesY[namsY],"[[","name")),collapse=","),")"),sep=" _ ")) 
    }))
  PostJoint <- data.frame(aggregate(JointModP$logpost,by=list(JointModP$ModNam),FUN = mean))
  names(PostJoint) <- c("Model","logpost")
  PostJoint$postP <- round(exp(PostJoint$logpost-max(PostJoint$logpost))/
                             sum(exp(PostJoint$logpost-max(PostJoint$logpost))),5)
  PostJoint$postF <- round(as.numeric(table(JointModP$ModNam))/nrow(JointModP),5)
  PostJoint <- PostJoint[order(PostJoint$postP,decreasing=T),]
  
  #-----
  
  
  nodesZ <- with(objVarSel,nodesZ)
  modnodeslZ = lapply(objVarSel$ModelStorageZ$model,"[[","nodes")
  names(modnodeslZ)=NULL
  MPIPsZp <- round(with(objVarSel,PostIncProb(modnodeslZ[orderz],postZp$post,nodesZ)),4)
  nodesY <- with(objVarSel,nodesY)
  modnodeslY = lapply(objVarSel$ModelStorageY$model,"[[","nodes")
  names(modnodeslY)=NULL
  MPIPsYp <- round(with(objVarSel,PostIncProb(modnodeslY[ordery],postYp$post,nodesY)),4)
  
  modnodeslZf = lapply(objVarSel$ModelStorageZ$model[names(postZf)],"[[","nodes")
  modnodeslYf = lapply(objVarSel$ModelStorageY$model[names(postYf)],"[[","nodes")
  names(modnodeslZf)=names(modnodeslYf)=NULL
  MPIPsZf <- round(with(objVarSel,PostIncProb(modnodeslZf,postZf,nodesZ)),4)
  MPIPsYf <- round(with(objVarSel,PostIncProb(modnodeslYf,postYf,nodesY)),4)
  
  postZf <- data.frame(models.z=with(objVarSel,unlist(lapply(names(postZf),function(x){
    nams=as.numeric(strsplit(x,split=",")[[1]])
    paste0(unlist(lapply(nodesZ[nams],"[[","name")),collapse=",")}))),post=postZf)[orderzf,c(1,3)]
  
  postYf <- data.frame(models.y=with(objVarSel,unlist(lapply(names(postYf),function(x){
    nams=as.numeric(strsplit(x,split=",")[[1]])
    paste0(unlist(lapply(nodesY[nams],"[[","name")),collapse=",")}))),post=postYf)[orderyf,c(1,3)]
  
  return(list(postYf=postYf,postZf=postZf,
              postYp=postYp,postZp=postZp,
              PostFJoint=PostFJoint,
              PostJoint=PostJoint,
              MPIPsYp=MPIPsYp,MPIPsZp=MPIPsZp,
              MPIPsYf=MPIPsYf,MPIPsZf=MPIPsZf,
              modnodeslZf=modnodeslZf))
}

#load("OutputBlueHawker20K.RData")
# ResHIP=summary.fun(testrun,priortype=NULL,erg.fn=mean)
# ResHOP=summary.fun(testrun,priortype="HOP")
# ResEPP=summary.fun(testrun,priortype="EPP")
# require(xtable)
# xtable(ResHIP$postZf)
# xtable(data.frame(pp=ResHIP$MPIPsYp,pf=ResHIP$MPIPsYf))
# xtable(data.frame(pp=ResHIP$MPIPsZp,pf=ResHIP$MPIPsZf))
# 
# 
# ResSimsHIP=summary.fun(simtestrun,priortype=NULL)
# ResSimsHOP=summary.fun(simtestrun,priortype="HOP")
# ResSimsEPP=summary.fun(simtestrun,priortype="EPP")
# require(xtable)
# xtable(ResSimsHIP$postZp)
# xtable(data.frame(pp=ResHIP$MPIPsYp,pf=ResSimsHIP$MPIPsYf))
# xtable(data.frame(pp=ResHIP$MPIPsZp,pf=ResHIP$MPIPsZf))
