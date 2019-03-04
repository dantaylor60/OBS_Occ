setwd(paste0(getwd(),"/PackageFunctions"))
source("OccupancySel_SH.R")

library(unmarked)
library(combinat)

data(mallard)
mallard.y[mallard.y>1]<-1
all0s <- which(apply(mallard.y,1,function(rr)sum(is.na(rr))==3))
Y <- mallard.y[-all0s,]

dataZ <- mallard.site[-all0s,]
mallard.obslist <- mallard.obs
mallard.obslist <- lapply(mallard.obslist,function(xx){
  xx = xx[-all0s,]
  return(xx)
  })
  
dataY <- data.frame(ivel = c(t(mallard.obslist$ivel)), 
                             date =  c(t(mallard.obslist$date)))
not.NA.obsY <- (1:(prod(dim(Y))))[!is.na(c(t(Y)))]
dataY <-dataY[not.NA.obsY,]


ff.Z<-as.formula(. ~ elev+length+forest+elev*length+elev*forest+elev^2 )
ff.Y<-as.formula(. ~ ivel + date + ivel^2 + date^2 + ivel*date)
start.formulaY<-as.formula(.~1)
start.formulaZ<-as.formula(.~1)
niter=10000

OBayesMallard <- Occup.IPSel(Y=Y,full.formulaZ=ff.Z,dataZ=dataZ,base.formulaZ=as.formula(.~ 1),full.formulaY=ff.Y,
                             dataY=dataY,base.formulaY=as.formula(.~ 1),SH=TRUE,maxdegZ=2,maxdegY=2,nodes.to.removeZ=NULL,nodes.to.removeY=NULL,
                             model.prior.type=c("HIP","EPP","HOP","HUP"),pxda=T,model.prior.pars="children",beta.prior.type="IP",
                             beta.prior.pars=list(alpha=1,nu=1),eff=TRUE,centerscale=TRUE,start.formulaZ=start.formulaZ,
                             start.formulaY=start.formulaY,niter=niter,ordy.by.site=T,search.start=T,
                             pshuffle=0.1,burnin=1000,thin=10)

source("SummaryCode.R")
ResHIP=summary.fun(OBayesMallard,priortype=NULL,erg.fn=median)
ResHOP=summary.fun(OBayesMallard,priortype="HOP")
ResHUP=summary.fun(OBayesMallard,priortype="HUP")
ResEPP=summary.fun(OBayesMallard,priortype="EPP")


ResHIP$MPIPsYf
ResHIP$MPIPsYp

ResHIP$MPIPsZf
ResHIP$MPIPsZp


indmodsY_p.in.f<-match(ResHIP$postYp[,1],ResHIP$postYf[,1])
indmodsZ_p.in.f<-match(ResHIP$postZp[,1],ResHIP$postZf[,1])

require(xtable)
xtable(data.frame(ResHIP$postZp,post.freq=ResHIP$postZf[indmodsZ_p.in.f,2]))
xtable(data.frame(ResHIP$postYp,post.freq=ResHIP$postYf[indmodsY_p.in.f,2]))


