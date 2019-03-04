psr<-function(full.formula,data,base.formula=as.formula(.~ 1),SH=TRUE,maxdeg=0,nodes.to.remove=NULL,model.prior.type="EPP",
              model.prior.pars="children",beta.prior.type="IP",beta.prior.pars=list(alpha=1,nu=1),eff=TRUE,centerscale=FALSE,
              start.formula=as.formula(.~ 1),niter=10000,ModelStorage,model.sequence,logpostvec,
              MatParRB,latrespini,shuffle=0){
  #full.formula is either a formula w/o parentheses
  #data is a data frame
  #base.formula is a formula
  #SH is a logical flag for strong heredity
  #maxdeg - if zero, then just interprets formula, if >0, then takes all vars in formula and makes surface of max degree in those vars
  #noderemoval vector of strings to remove (we will still force WFM on the full model)
  #msp, mspp, pp, ppp, eff are the usual
  #centerscale - center and scale the vars - defaults to false
  #starting point - formula - always becomes WF - if not given, defaults to base model
  #niter - numer of iterations for the stochastic search algorithm
  #respname - name of the response as coded inside of the data.frame "data"
 
  
  #-------------------------------------------------------------------------
{
{

  if(centerscale){
    isfact = names(unlist(lapply(data,class)))[(unlist(lapply(data,class)))!="factor"]
    data[,isfact] = scale(data[,isfact])
  }
  full.get<-function(formula,intercept.label="1"){
    vars<-all.vars(formula)
    vars<-vars[2:length(vars)]
    Lv<-length(vars)
    nodes<-as.character(formula)[3]
    nodes<-gsub(pattern="\n", replacement="",x=nodes)
    nodes<-gsub(pattern=" ", replacement="",x=nodes)
    if(grepl("(",nodes,fixed=TRUE)==TRUE){return("Formula must be formed without parentheses.")}
    nodes<-strsplit(nodes,"+",fixed=TRUE)[[1]]
    Ln<-length(nodes)
    out<-matrix(0,Ln,Lv)
    colnames(out)<-vars
    for(i in 1:Ln){
      node.loc<-strsplit(nodes[i],"*",fixed=TRUE)[[1]]
      Ln.loc<-length(node.loc)
      for(j in 1:Ln.loc){
        var.loc<-strsplit(node.loc[j],"^",fixed=TRUE)[[1]]
        if(length(var.loc)>1){
          out[i,var.loc[1]]<-out[i,var.loc[1]]+as.integer(var.loc[2])
        }else{
          out[i,var.loc[1]]<-out[i,var.loc[1]]+1
        }
      }
    }
    out<-unique(out,MARGIN=1)
    i<-1
    while(i<=dim(out)[1]){
      ind<-which(out[i,]>0)
      Li<-length(ind)
      if(Li>0){
        sub<-matrix(out[i,],Li,Lv,byrow=TRUE)
        for(j in 1:Li){
          sub[j,ind[j]]<-sub[j,ind[j]]-1
        }
        out<-rbind(out,sub)
      }
      i<-i+1
    }
    out<-unique(out,MARGIN=1)
    degree<-apply(out,1,sum)
    degree<-sort(degree,index.return=TRUE)
    
    out<-matrix(out[degree$ix,],ncol=Lv)
    colnames(out)<-vars
    degree<-degree$x
    namer<-function(x,labels,intercept.label){
      ind<-which(x>0)
      Li<-length(ind)
      if(Li==0){return(as.character(intercept.label))}else{
        if(Li==1){
          if(x[ind[1]]==1){
            return(labels[ind[1]])
          }else{
            return(paste(labels[ind[1]],"^",x[ind[1]],sep=""))
          }
        }else{
          out<-""
          for(i in 1:(Li-1)){
            if(x[ind[i]]==1){
              out<-paste(out,labels[ind[i]],"*",sep="")
            }else{
              out<-paste(out,labels[ind[i]],"^",x[ind[i]],"*",sep="")
            }  			
          }
          if(x[ind[Li]]==1){
            out<-paste(out,labels[ind[Li]],sep="")
          }else{
            out<-paste(out,labels[ind[Li]],"^",x[ind[Li]],sep="")
          }				
          return(out)
        }
      }
    }
    
    out.vec<-apply(out,1,namer,labels=colnames(out),intercept.label=intercept.label)
    
    rownames(out)<-out.vec
    
    get.edges<-function(out,degree){
      edges<-cbind(c(rownames(out)[1]),rownames(out)[degree==1])
      if(max(degree)==1){
        return(edges)
      }else{
        comp<-function(x,y){
          sum(abs(x-y))
        }
        for(j in 1:(max(degree)-1)){
          vertices.loc<-rownames(out)[degree==j]
          d<-j+1
          ind.loc<-which(degree==d)
          out.loc<-matrix(out[ind.loc,],ncol=dim(out)[2])
          rownames(out.loc)<-rownames(out)[ind.loc]
          for(vertex in vertices.loc){
            points.to<-apply(out.loc,1,comp,out[vertex,])
            points.to<-which(points.to==1)
            if(length(points.to)>0){
              edges<-rbind(edges,cbind(c(vertex),rownames(out.loc)[points.to]))
            }
          }
        }
        return(edges)
      }
    }
    edges<-get.edges(out,degree)
    
    get.parents<-function(out.vec,edges){
      parents<-list()
      for(i in 1:length(out.vec)){
        parents[[out.vec[i]]]<-edges[edges[,2]==out.vec[i],1]
      }
      parents
    }
    parents<-get.parents(out.vec,edges)
    
    return(list(vertices=rownames(out),nodes=out,edges=edges,degree=degree,parents=parents))
  }
  if(maxdeg==0){
    MF.mat<-full.get(full.formula)$nodes
    vars<-colnames(MF.mat)
    colnames(MF.mat)<-NULL
  }else{
    vars<-all.vars(full.formula)
    vars<-vars[2:length(vars)]
    
    require(combinat)
    allterm.generator <-function(v,r){
      do.call(rbind,lapply(1:r,FUN=function(v,x){
        t(xsimplex(p=v,n=x))},v=v))
    }
    
    MF.mat<-rbind(c(0),allterm.generator(length(vars),maxdeg))
    namer<-function(x,labels,intercept.label){
      ind<-which(x>0)
      Li<-length(ind)
      if(Li==0){return(as.character(intercept.label))}else{
        if(Li==1){
          if(x[ind[1]]==1){
            return(labels[ind[1]])
          }else{
            return(paste(labels[ind[1]],"^",x[ind[1]],sep=""))
          }
        }else{
          out<-""
          for(i in 1:(Li-1)){
            if(x[ind[i]]==1){
              out<-paste(out,labels[ind[i]],"*",sep="")
            }else{
              out<-paste(out,labels[ind[i]],"^",x[ind[i]],"*",sep="")
            }    		
          }
          if(x[ind[Li]]==1){
            out<-paste(out,labels[ind[Li]],sep="")
          }else{
            out<-paste(out,labels[ind[Li]],"^",x[ind[Li]],sep="")
          }				
          return(out)
        }
      }
    }
    rownames(MF.mat)<-apply(MF.mat,1,namer,labels=vars,intercept.label="1")
    
    
  }
}
{
  base.get<-function(formula,intercept.label="1"){
    vars<-all.vars(formula)
    if(length(vars)==1){
      return(intercept.label)
    }else{
      return(full.get(formula,intercept.label)$vertices)
    }
  }
}
{
  base.nums<-base.get(base.formula)
  MB.mat<-matrix(MF.mat[base.nums,],ncol=dim(MF.mat)[2])
  start.nums<-base.get(start.formula)
  M.start.mat<-matrix(MF.mat[start.nums,],ncol=dim(MF.mat)[2])
  
  node.names<-row.names(MF.mat)
  rownames(MF.mat)<-NULL
  
  ## frame for X matrix
  X.frame<-data.frame(data[,vars],row.names=row.names(data))
  names(X.frame)=vars
  vars.classes<-sapply(X.frame,class)
  fac.ind<-vars.classes=="factor"
  names(vars.classes)<-NULL
  
  if(sum(fac.ind)>1){
    to.remove<-which(apply(MF.mat,1,function(x){max(x[fac.ind])>1}))
    MF.mat<-MF.mat[-to.remove,]
    node.names<-node.names[-to.remove]
    if(dim(MB.mat)[1]>1){
      to.remove<-which(apply(MB.mat,1,function(x){max(x[fac.ind])>1}))
      MB.mat<-MB.mat[-to.remove,]
    }
    to.remove<-which(apply(M.start.mat,1,function(x){max(x[fac.ind])>1}))
    M.start.mat<-M.start.mat[-to.remove,]
  }
}  
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
## Given MF and MB matrices
## we set up the the matrices to coincide with base nodes being at the top
{  
  MF.deg<-apply(MF.mat,1,sum)
  MF.sortrows<-sort(MF.deg,decreasing=FALSE,index.return=TRUE)
  MF.mat<-MF.mat[MF.sortrows$ix,]
  MF.deg<-MF.sortrows$x
  node.names<-node.names[MF.sortrows$ix]
  MB.find<-function(b){
    a<-apply(MF.mat,1,identical,y=b)==1
    if(prod(!a)){
      0
    }else{
      which(a)
    }
  }
  base.nodes<-sort(apply(MB.mat,1,MB.find))
  if(min(base.nodes)==0){
    return("Base Model is not contained in Full Model.")
  }
  MB.mat<-matrix(MF.mat[base.nodes,],ncol=dim(MF.mat)[2])
  rownames(MB.mat)<-rownames(MF.mat)[base.nodes]
  MF.mat<-rbind(MF.mat[base.nodes,],MF.mat[-base.nodes,])
  rownames(MF.mat)[1:length(base.nodes)]<-rownames(MB.mat)
  MF.deg<-c(MF.deg[base.nodes],MF.deg[-base.nodes])
  L.base.nodes<-length(base.nodes)
  node.names<-c(node.names[base.nodes],node.names[-base.nodes])
  base.nodes<-c(1:L.base.nodes)
}
## MF.mat, MB.mat, base.nodes, and L.base.nodes are assumed to exist from now on
## as is MF.deg
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
## We now build the nodes list
## Each node is a list of 
## vec from MF.mat
## parents
## children
## degree
## name
## formula
#-------------------------------------------------------------------------
{
  DIM<-dim(MF.mat)
  rn<-c(1:DIM[1])
  rownames(MF.mat)<-rn
  rownames(MB.mat)<-c(1:L.base.nodes)
  cn<-c(1:DIM[2])
  nodes<-vector(mode="list",length=DIM[1])
  for(i in rn){
    nodes[[i]]<-list(vec=MF.mat[i,],children=NULL,parents=NULL)
  }
  Eye<-diag(DIM[2])
  v.c.eq<-function(v,c=Eye,d=1){
    as.integer(rownames(c)[which(apply(c,d,identical,y=v)==TRUE)])
  }
  for(i in MF.deg[1]:(MF.deg[DIM[1]]-1)){
    lo<-which(MF.deg==i)
    hi<-which(MF.deg==(i+1))
    HM<-matrix(MF.mat[hi,],ncol=dim(MF.mat)[2])
    rownames(HM)<-rownames(MF.mat)[hi]
    for(j in lo){
      V<-MF.mat[j,]+Eye
      nctemp<-unlist(apply(V,2,v.c.eq,c=HM))
      if(length(nctemp)>0){
        nodes[[j]]$children<-nctemp
        for(k in nodes[[j]]$children){
          nodes[[k]]$parents<-c(nodes[[k]]$parents,j)
        }
      }
    }
  }
  make.formula<-function(vec){
    ind<-(vec>0)
    if(sum(ind)==0){
      return("1")
    }
    V<-vars[ind]  
    C<-vars.classes[ind]
    vec<-vec[ind]
    ind<-(C=="factor")
    
    if(prod(ind)){
      out<-paste(V,collapse=":")
    }else if(prod(!ind)){
      out<-paste("I(",paste(paste(V,vec,sep="^"),collapse="*"),")",sep="")
    }else{
      a<-paste("I(",paste(paste(V[!ind],vec[!ind],sep="^"),collapse="*"),")",sep="")
      b<-paste(V[ind],collapse=":")
      out<-paste(a,b,sep=":")
    }
    out
  }
  
  
  for(i in rn){
    nodes[[i]]$degree<-as.numeric(MF.deg[i])
    nodes[[i]]$name<-node.names[i]
    nodes[[i]]$formula<-make.formula(nodes[[i]]$vec)
  }
  L.nodes<-length(nodes)
}  
#-------------------------------------------------------------------------
## "nodes" is assumed to exist from here on out, as is L.nodes


#-------------------------------------------------------------------------
## Now we consider node removal
{  
  if(!is.null(nodes.to.remove)){
    fREM<-formula(paste(".~",paste(nodes.to.remove,collapse="+"),collapse=""))
    remover<-function(formula){
      Lv<-length(vars)
      nodes<-as.character(formula)[3]
      if(grepl("(",nodes,fixed=TRUE)==TRUE){return("Formula must be formed without parentheses.")}
      nodes<-strsplit(nodes," + ",fixed=TRUE)[[1]]
      Ln<-length(nodes)
      out<-matrix(0,Ln,Lv)
      colnames(out)<-vars
      for(i in 1:Ln){
        node.loc<-strsplit(nodes[i]," * ",fixed=TRUE)[[1]]
        Ln.loc<-length(node.loc)
        for(j in 1:Ln.loc){
          var.loc<-strsplit(node.loc[j],"^",fixed=TRUE)[[1]]
          if(length(var.loc)>1){
            out[i,var.loc[1]]<-out[i,var.loc[1]]+as.integer(var.loc[2])
          }else{
            out[i,var.loc[1]]<-out[i,var.loc[1]]+1
          }
        }
      }
      out<-unique(out,MARGIN=1)
      out
    }
    M.rem.mat<-remover(fREM)
    colnames(M.rem.mat)<-NULL
    rem.nodes<-apply(M.rem.mat,1,function(x){which(apply(MF.mat,1,identical,y=x))})
    rem.nodes<-rem.nodes[!(rem.nodes%in%base.nodes)]
    rem.nodes<-list(rem.nodes)
    i<-2
    flag.rem<-1
    while(flag.rem){
      rem.nodes.new<-unique(unlist(lapply(nodes[rem.nodes[[i-1]]],function(n){n$children})))
      flag.rem<-length(rem.nodes.new)>0
      if(flag.rem){
        rem.nodes[[i]]<-rem.nodes.new
        i<-i+1 
      }
    }
    rem.nodes<-unlist(rem.nodes)
    nodes.map<-rep(NA,L.nodes)
    nodes.map[-rem.nodes]<-c(1:(L.nodes-length(rem.nodes)))
    MF.mat<-MF.mat[-rem.nodes,]
    nodes<-nodes[-rem.nodes]
    MF.deg<-MF.deg[-rem.nodes]
    node.names<-node.names[-rem.nodes]
    L.nodes<-length(nodes)
    rownames(MF.mat)<-c(1:L.nodes)
    for(i in 1:L.nodes){
      nodes[[i]]$children<-nodes.map[nodes[[i]]$children]
      nodes[[i]]$children<-nodes[[i]]$children[!is.na(nodes[[i]]$children)]
      if(length(nodes[[i]]$children)==0){
        nodes[[i]]<-list(vec=nodes[[i]]$vec,children=NULL,parents=nodes[[i]]$parents,degree=nodes[[i]]$degree,name=nodes[[i]]$name)
      }
    }
  }
}
#-------------------------------------------------------------------------


#------------------------------------------------------------------------- 
## We have to have a beginning model and have to do four things to it
## First, we have to check if it is the base or full model
## Second, we have to check if it satisfies our heredity condition
## Third, we have to find its extreme nodes
## Fourth, we have to find its children nodes
## The starting model is called M in these functions and assumed to only have M$nodes
## It is called M.start in the broader code
{
  M.start<-list(nodes=sort(unlist(apply(M.start.mat,1,function(x){which(apply(MF.mat,1,identical,y=x))}))))
  names(M.start$nodes)<-NULL
  
  Model.check.base<-function(M){
    as.logical(prod(base.nodes%in%M$nodes))
  }
  
  Model.check.heredity<-function(M,SH=TRUE){
    if(!(identical(base.nodes,M$nodes))){
      if(SH){
        as.logical(prod(unlist(lapply(nodes[M$nodes[-base.nodes]],function(n){prod(n$parents%in%M$nodes)==1}))))
      }else{
        as.logical(prod(unlist(lapply(nodes[M$nodes[-base.nodes]],function(n){sum(n$parents%in%M$nodes)>0}))))
      }
    }else{
      TRUE
    }
  }
  
  Model.get.extreme<-function(M,SH=TRUE){
    if(SH){
      out<-M$nodes[-base.nodes][unlist(lapply(nodes[M$nodes[-base.nodes]],function(n){sum(n$children%in%M$nodes)==0}))]
    }else{
      cnodes<-lapply(nodes[M$nodes[-base.nodes]],function(n){n$children[n$children%in%M$nodes]})
      to.add.logical<-unlist(lapply(cnodes,function(x){length(x)==0}))
      ind<-which(!to.add.logical)
      if(!prod(to.add.logical)){
        to.add.logical[ind]<-unlist(lapply(cnodes[ind],function(cn){as.logical(prod(unlist(lapply(nodes[cn],function(n){sum(n$parents%in%M$nodes)>1}))))}))
      }
      out<-M$nodes[-base.nodes][to.add.logical]
    }
    if(length(out)==0){NULL}else{out}
  }
  
  Model.get.children<-function(M,SH=TRUE){
    out<-unique(sort(unlist(lapply(nodes[M$nodes],function(n){n$children[!(n$children%in%M$nodes)]}))))
    if(SH){
      out<-out[unlist(lapply(nodes[out],function(n){as.logical(prod(n$parents%in%M$nodes))}))]
    }
    out
  }
  ## Now we get the children of the base model and the extereme of the full model
  MB.children.fun<-function(SH=TRUE){
    if(SH){
      out<-which(unlist(lapply(nodes,function(n){prod(n$parents%in%base.nodes)==1})))
    }else{
      out<-which(unlist(lapply(nodes,function(n){sum(n$parents%in%base.nodes)>0})))
    }
    out<-out[out>L.base.nodes]
    out
  }
  
  MF.extreme.fun<-function(SH=TRUE){
    if(SH){
      out<-which(as.vector(lapply(nodes,function(l){is.null(l$children)}))==TRUE)
    }else{
      out<-which(as.vector(lapply(nodes,function(l){length(l$parents)==1}))==TRUE)
      out<-unlist(lapply(nodes[out],function(l)(l$parents)))
      out<-c(1:length(nodes))[-out]
    }
    out<-out[out>L.base.nodes]
    out
  }
}

#------------------------------------------------------------------------- 

#-------------------------------------------------------------------------
#---- Define mutable object for model storage
{
  ## We start out ModelStorage$model by including both the base and full model
  ModelStorage$model[[1]]<-list(nodes=base.nodes,extreme=NULL,children=MB.children.fun(SH))
  names(ModelStorage$model)[1]<-paste(base.nodes,collapse=",")
  ModelStorage$model[[paste(c(1:L.nodes),collapse=",")]]<-list(nodes=c(1:L.nodes),extreme=MF.extreme.fun(SH),children=NULL)
}
#-------------------------------------------------------------------------
## This is our model store, its names are given by
## names(ModelStorage$model)[i]=paste(models[[i]]$nodes,collapse=",")
## This is now assumed to exist and will be dynamic

#-------------------------------------------------------------------------
#Base, Full and Start model checks
{
  ## Check Heredity for Base Model
  if(!Model.check.heredity(ModelStorage$model[[1]],SH)){
    return("Base Model does not satisfy Heredity Condition.")
  }else{
    ModelStorage$model[[1]]$children<-MB.children.fun(SH)
  }
  ## Check Heredity for Full Model
  if(!Model.check.heredity(ModelStorage$model[[2]],SH)){
    return("Full Model does not satisfy Heredity Condition.")
  }else{
    ModelStorage$model[[2]]$extreme<-MF.extreme.fun(SH)
  }
  
  ## Make Sure that Base is in Start
  if(!Model.check.base(M.start)){
    M.start$nodes<-sort(unique(c(M.start$nodes,base.nodes)))
  }
  ## If Start is Base or Full, just go to it
  if(identical(M.start$nodes,base.nodes)){
    M.start<-ModelStorage$model[[1]]
  }else if(identical(M.start$nodes,ModelStorage$model[[2]]$nodes)){
    M.start<-ModelStorage$model[[2]]
  }else{
    if(!Model.check.heredity(M.start,SH)){
      M.start<-ModelStorage$model[[1]]
    }else{
      M.start$extreme<-Model.get.extreme(M.start,SH)
      M.start$children<-Model.get.children(M.start,SH)
      ModelStorage$model[[3]]<-M.start
      names(ModelStorage$model)[3]<-paste(ModelStorage$model[[3]]$nodes,collapse=",")
    }
  }
}
#-------------------------------------------------------------------------
##Now we have a starting model and models  


#-------------------------------------------------------------------------
#define functions used to calculate logPriors  
{
  ## For the model space prior
  modelPriormaker <- function(MF.deg,type,hyperpars=c("unif","children"),base.nodes){
    if(type=="EPP"){logprior<-function(M){0}}
    alpha <- 1
    L.base.nodes<-length(base.nodes)
    L.MF<-length(MF.deg)
    if(hyperpars=="unif"){ ## No children
      beta <- 1
      lbab<-lbeta(alpha,beta)
      if(type=="HUP"){
        logprior<-function(M){
          lbeta(length(M$nodes)-L.base.nodes+alpha,length(M$children)+beta)-lbab
        }
      }else if(type=="HOP"){
        logprior<-function(M){
          nodesnotbase <- M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          G<-table(factor(deg.nnb,levels=c(mindeg:maxdeg)))
          CH<-table(factor(deg.ch,levels=c(mindeg:maxdeg)))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          sum(lbeta(G+alpha,CH+beta)-lbab)
        }
      }else if (type=="HIP"){
        lprob<-log(alpha)-log(alpha+beta)
        lprobm<-log(beta)-log(alpha+beta)
        logprior<-function(M){
          nodesnotbase <- M$nodes[-base.nodes]
          length(nodesnotbase)*lprob+length(M$children)*lprobm
        }
      }else if (type=="HLP"){
        logprior<-function(M){
          nodesnotbase<-M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          pfun<-function(n){sum(n$vec>0)}
          p.nnb<-unlist(lapply(nodes[nodesnotbase],pfun))
          p.ch<-unlist(lapply(nodes[M$children],pfun))
          minp<-min(c(p.nnb,p.ch))
          maxp<-max(c(p.nnb,p.ch))
          G<-as.vector(xtabs(~factor(deg.nnb,levels=c(mindeg:maxdeg))+factor(p.nnb,levels=c(minp:maxp))))
          CH<-as.vector(xtabs(~factor(deg.ch,levels=c(mindeg:maxdeg))+factor(p.ch,levels=c(minp:maxp))))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          sum(lbeta(G+alpha,CH+beta)-lbab)
        }
      }else if (type=="HPP"){
        logprior<-function(M){
          nodesnotbase<-M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          pfun<-function(n){sum(!(n$parents%in%base.nodes))}
          p.nnb<-unlist(lapply(nodes[nodesnotbase],pfun))
          p.ch<-unlist(lapply(nodes[M$children],pfun))
          minp<-min(c(p.nnb,p.ch))
          maxp<-max(c(p.nnb,p.ch))
          G<-as.vector(xtabs(~factor(deg.nnb,levels=c(mindeg:maxdeg))+factor(p.nnb,levels=c(minp:maxp))))
          CH<-as.vector(xtabs(~factor(deg.ch,levels=c(mindeg:maxdeg))+factor(p.ch,levels=c(minp:maxp))))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          sum(lbeta(G+alpha,CH+beta)-lbab)
        }
      }else if (type=="HTP"){
        logprior<-function(M){
          nodesnotbase<-M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          pfun<-function(n){paste(sort(n$vec[n$vec>0]),collapse="'")}
          p.nnb<-unlist(lapply(nodes[nodesnotbase],pfun))
          p.ch<-unlist(lapply(nodes[M$children],pfun))
          lev.p<-unique(c(unique(p.nnb),unique(p.ch)))
          G<-as.vector(xtabs(~factor(deg.nnb,levels=c(mindeg:maxdeg))+factor(p.nnb,levels=lev.p)))
          CH<-as.vector(xtabs(~factor(deg.ch,levels=c(mindeg:maxdeg))+factor(p.ch,levels=lev.p)))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          sum(lbeta(G+alpha,CH+beta)-lbab)
        }
      }
    }else{ ## here we use children
      if(type=="HUP"){
        beta <- L.MF-L.base.nodes
        lbab<-lbeta(alpha,beta)
        logprior<-function(M){
          lbeta(length(M$nodes)-L.base.nodes+alpha,length(M$children)+beta)-lbab
        }
      }else if(type=="HOP"){
        logprior<-function(M){
          nodesnotbase <- M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          G<-table(factor(deg.nnb,levels=c(mindeg:maxdeg)))
          CH<-table(factor(deg.ch,levels=c(mindeg:maxdeg)))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          beta<-G+CH
          sum(lbeta(G+alpha,CH+beta)-lbeta(alpha,beta))
        }
      }else if(type=="HIP"){
        logprior<-function(M){
          nodesnotbase <- M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          G<-table(factor(deg.nnb,levels=c(mindeg:maxdeg)))
          CH<-table(factor(deg.ch,levels=c(mindeg:maxdeg)))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          beta<-G+CH
          lbab<-log(alpha+beta)
          sum(G*(log(alpha)-lbab)+CH*(log(beta)-lbab))
        }
      }else if (type=="HLP"){
        logprior<-function(M){
          nodesnotbase<-M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          pfun<-function(n){sum(n$vec>0)}
          p.nnb<-unlist(lapply(nodes[nodesnotbase],pfun))
          p.ch<-unlist(lapply(nodes[M$children],pfun))
          minp<-min(c(p.nnb,p.ch))
          maxp<-max(c(p.nnb,p.ch))
          G<-as.vector(xtabs(~factor(deg.nnb,levels=c(mindeg:maxdeg))+factor(p.nnb,levels=c(minp:maxp))))
          CH<-as.vector(xtabs(~factor(deg.ch,levels=c(mindeg:maxdeg))+factor(p.ch,levels=c(minp:maxp))))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          beta<-G+CH
          sum(lbeta(G+alpha,CH+beta)-lbeta(alpha,beta))
        }
      }else if (type=="HPP"){
        logprior<-function(M){
          nodesnotbase<-M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          pfun<-function(n){sum(!(n$parents%in%base.nodes))}
          p.nnb<-unlist(lapply(nodes[nodesnotbase],pfun))
          p.ch<-unlist(lapply(nodes[M$children],pfun))
          minp<-min(c(p.nnb,p.ch))
          maxp<-max(c(p.nnb,p.ch))
          G<-as.vector(xtabs(~factor(deg.nnb,levels=c(mindeg:maxdeg))+factor(p.nnb,levels=c(minp:maxp))))
          CH<-as.vector(xtabs(~factor(deg.ch,levels=c(mindeg:maxdeg))+factor(p.ch,levels=c(minp:maxp))))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          beta<-G+CH
          sum(lbeta(G+alpha,CH+beta)-lbeta(alpha,beta))
        }
      }else if (type=="HTP"){
        logprior<-function(M){
          nodesnotbase<-M$nodes[-base.nodes]
          deg.nnb<-MF.deg[M$nodes][-base.nodes]
          deg.ch<-MF.deg[M$children]
          mindeg<-min(c(deg.nnb,deg.ch))
          maxdeg<-max(c(deg.nnb,deg.ch))
          pfun<-function(n){paste(sort(n$vec[n$vec>0]),collapse="'")}
          p.nnb<-unlist(lapply(nodes[nodesnotbase],pfun))
          p.ch<-unlist(lapply(nodes[M$children],pfun))
          lev.p<-unique(c(unique(p.nnb),unique(p.ch)))
          G<-as.vector(xtabs(~factor(deg.nnb,levels=c(mindeg:maxdeg))+factor(p.nnb,levels=lev.p)))
          CH<-as.vector(xtabs(~factor(deg.ch,levels=c(mindeg:maxdeg))+factor(p.ch,levels=lev.p)))
          Z<-(G==0)*(CH==0)
          G<-G[!Z]
          CH<-CH[!Z]
          beta<-G+CH
          sum(lbeta(G+alpha,CH+beta)-lbeta(alpha,beta))
        }
      }
    }
    return(logprior)
  }
  modelPriormakerV <- Vectorize(FUN=modelPriormaker,vectorize.args="type")
  modpriorfn<-function(M){
    unlist(lapply(modelPriormakerV(MF.deg=MF.deg,type=model.prior.type,hyperpars=model.prior.pars,base.nodes=base.nodes),
           function(prfn)prfn(M)))
  }
}
#---- Calculates model fit and generate MatPar
{
  baseModinfo<-function(mat0,mat){
    n<-dim(mat)[1]
    mat0<-matrix(mat0,nrow=n)
    pbase<-dim(mat0)[2]
    namesFull = colnames(mat)
    pmx<-length(namesFull)
    xmat0inv<-armaSolve(crossprod(mat0))
    Hmat0 <-mat0%*%tcrossprod(xmat0inv,mat0)
    return(list(mat0=mat0,namesFull=namesFull,Hmat0=Hmat0,n=n,pbase=pbase,pmx=pmx,
                xmat0inv=xmat0inv))
  }
  formBase = unlist(lapply(nodes[ModelStorage$model[[1]]$nodes],"[","formula"))
  formFull = unlist(lapply(nodes[ModelStorage$model[[2]]$nodes],"[","formula"))
  names(formBase)=names(formFull)=NULL
  formFull=formFull[!(formFull%in%formBase)]
  matM0 <-  model.matrix(terms(reformulate(formBase),keep.order = TRUE),X.frame)
  matMF <- model.matrix(terms(reformulate(formFull),
                              keep.order = TRUE),X.frame)[,-1]
  baseinfo <- baseModinfo(mat0=matM0,mat=matMF)
  
  fn.fit<- function(M,latresp){
    form.names<- unlist(lapply(nodes[M$nodes],"[","formula"))
    names(form.names)=NULL
    out<-with(baseinfo,{
      #--- pars base model
      latresp = matrix(latresp,ncol=1)
      parRB <- params <-rep(0,(pmx+pbase))
      parRB[1:pbase] <- xmat0inv%*%t(mat0)%*%latresp
      params[1:pbase] <- rmvnormArma(1, parRB[1:pbase], xmat0inv)
      bcst <- as.numeric(crossprod(latresp,mat0)%*%parRB[1:pbase])
      if(!identical(form.names,formBase)){
        form.model=form.names[!(form.names%in%formBase)]
        matCand <- model.matrix(terms(reformulate(form.model),keep.order = TRUE),X.frame)
        namMatCand=dimnames(matCand)[[2]][-1]
        matCand = as.matrix(matCand[,-1])
        colnames(matCand)= namMatCand
        
        CoefCandinMF <- (pbase+which(namesFull%in%colnames(matCand)))#position cols MCand in MF
        sizeA = length(CoefCandinMF)
        pA = pbase+sizeA
        
        #--- pars candidate model
        matorth<-matCand-crossprod(Hmat0,matCand)
        xmatAo = crossprod(matCand)+(pA/(2*n))*crossprod(matorth)
        xmatAoinv = armaSolve(xmatAo)
        
        HAr = matorth%*%tcrossprod(xmatAoinv,matorth)
        
        parRB[CoefCandinMF] <- tcrossprod(xmatAoinv,matorth)%*%latresp
        
        params[CoefCandinMF] = rmvnormArma(1, parRB[CoefCandinMF], xmatAoinv)
        loglik = 0.5*(sizeA*log(pA/(2*n+pA))-crossprod(latresp,(diag(n)-Hmat0-HAr))%*%latresp)
        bcst <- as.numeric(bcst+crossprod(latresp,matCand)%*%parRB[CoefCandinMF])
        
      }else{
        loglik = -0.5*(crossprod(latresp,(diag(n)-Hmat0))%*%latresp)
        CoefCandinMF=NULL
        bcst <- as.numeric(crossprod(latresp,mat0)%*%parRB[c(1:pbase)])
      }
      return(list(parRB=parRB,params=params,loglik=loglik,bcst=bcst,
                  CoefCandinMF=c(1:pbase,CoefCandinMF) ) )
    })
    out
  }
  #----insert changes here with functions from CppFns.cpp
  ModFit.fn <- function(M,latresp){
    modelfitfn <- fn.fit(M,latresp)
    M$parsRB <- modelfitfn$parRB
    M$pars <- modelfitfn$params
    M$loglik <- modelfitfn$loglik
    M$modpr <-  unlist(modpriorfn(M))
    M$bcst <- modelfitfn$bcst
    return(M)
  }
  #Vectorize function
  VModFit.fn <- Vectorize(ModFit.fn,vectorize.args="M",SIMPLIFY=F)
  MatParRB$ParRB <- matrix(NA,nrow=niter,ncol=(ncol(matM0)+ncol(matMF)))
}
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
#calculate model fit for base, full and starting models
{
  BaseFitMod <- ModFit.fn(M=ModelStorage$model[[1]],latrespini)
  ModelStorage$model[[1]]<- BaseFitMod
  names(ModelStorage$model)[1] = paste(ModelStorage$model[[1]]$nodes,collapse=",")
  if(length(ModelStorage$model)==2){
    ModelStorage$model[[2]] <- ModFit.fn(ModelStorage$model[[2]],latresp=latrespini)
    names(ModelStorage$model)[2] = paste(ModelStorage$model[[2]]$nodes,collapse=",")
  }else if(length(ModelStorage$model)==3){
    ModelStorage$model[2:3] <- VModFit.fn(ModelStorage$model[2:3],latresp=latrespini)
    names(ModelStorage$model)[2:3] = unlist(lapply(ModelStorage$model[2:3],function(x)paste(x$nodes,collapse=",")))
  }
  numpost <- lapply(ModelStorage$model,function(x)(x$loglik+x$modpr[1]))
  names(numpost)=NULL
  logpostvec$logpost <- numpost
  Model.current<-ModelStorage$model[[length(ModelStorage$model)]]
}
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#Define functions used for the stochastic search algorithm
{  
  ## Functions for adding and subtracting nodes, also vectorized
  Model.add.node<-function(M,node,SH=TRUE){
    out<-M
    out$nodes<-sort(c(out$nodes,node))
    if(SH){
      out$extreme<-sort(c(out$extreme[!(out$extreme%in%nodes[[node]]$parents)],node))
      out$children<-out$children[out$children!=node]
      out$children<-sort(unique(
        c(out$children,
          nodes[[node]]$children[
            unlist(lapply(
              nodes[nodes[[node]]$children],
              function(l){as.logical(prod(l$parents%in%out$nodes))}
            )
            )
            ]
        )
      ))
    }else{
      pnode<-nodes[[node]]$parents
      pnode<-pnode[pnode%in%M$nodes]
      pnodeE<-pnode[pnode%in%M$extreme]
      if(length(pnode)==1 && length(pnodeE)==1){
        out$extreme<-out$extreme[-which(out$extreme==pnode)]
      }
      out$extreme<-sort(c(out$extreme,node))
      
      D<-nodes[[node]]$children
      D<-D[D%in%out$nodes]
      if(length(D)>0){
        Sigma<-unique(unlist(lapply(nodes[D],function(l){l$parents[l$parents%in%M$nodes]})))
        Sigma<-Sigma[!(Sigma%in%M$extreme)]
        if(length(Sigma)>0){
          Z<-lapply(nodes[Sigma],function(l){l$children[!(l$children%in%D) & (l$children%in%M$nodes)]})
          lZ<-lapply(Z,length)
          ind<-(lZ==0)
          Sigma.to.add.1<-Sigma[ind]
          if(length(Sigma)>0){
            Sigma<-Sigma[!ind]
            Z<-Z[!ind]
            ind<-unlist(lapply(Z,function(S){
              as.logical(prod(unlist(lapply(nodes[S],function(l){length(l$parents[l$parents%in%M$nodes])>1}))))
              }))
          
            Sigma.to.add.2<-Sigma[ind]  
            out$extreme<-c(out$extreme,Sigma.to.add.2)
          }
          out$extreme<-c(out$extreme,Sigma.to.add.1)
        }
      }
      out$children<-sort(unique(c(out$children[-which(out$children==node)],
                                  nodes[[node]]$children[!(nodes[[node]]$children%in%out$nodes)])))
    }
    out$extreme<-sort(unique(out$extreme[out$extreme>L.base.nodes]))
    out$children<-sort(unique(out$children[out$children>L.base.nodes]))
    out$loglik<-out$modpr<-out$parsRB<-out$pars<-NA
    out
  }
  V.Model.add.node<-Vectorize(Model.add.node,vectorize.args="node",SIMPLIFY=FALSE)
  
  Model.subtract.node<-function(M,node,SH=TRUE){
    out<-M
    out$nodes<-out$nodes[-which(out$nodes==node)]
    if(SH){
      out$extreme<-sort(
        c(out$extreme[-which(out$extreme==node)],
          nodes[[node]]$parents[
            unlist(lapply(
              nodes[nodes[[node]]$parents],
              function(l){as.logical(prod(!l$children%in%out$nodes))}
            ))
            ]
        )
      )
      out$children<-out$children[!(out$children%in%nodes[[node]]$children)]
      out$children<-sort(c(out$children,node))
    }else{
      C<-nodes[[node]]$children
      C<-C[C%in%M$children]
      ind<-which(unlist(lapply(nodes[C],function(l){
        length(l$parents[l$parents%in%out$nodes])
      }))==0)
      out$children<-sort(unique(c(out$children[!(out$children%in%C[ind])],node)))
      
      out$extreme<-out$extreme[-which(out$extreme==node)]
      
      pnode<-nodes[[node]]$parents
      pnode<-pnode[pnode%in%out$nodes]
      if(length(pnode)==1){
        cpnode<-nodes[[pnode]]$children
        cpnode<-cpnode[cpnode%in%out$nodes]
        pcpnode<-lapply(nodes[cpnode],function(l){l$parents[l$parents%in%out$nodes]})
        lpcpnode<-unlist(lapply(pcpnode,length))
        if(prod(lpcpnode>1)){
          out$extreme<-sort(c(out$extreme,pnode))
        }
      }
      
      C<-nodes[[node]]$children
      C<-C[C%in%M$nodes]
      B<-unlist(lapply(nodes[C],function(l){a<-l$parents[l$parents%in%out$nodes];a[1]*(length(a)==1)}))
      B<-B[B!=0]
      out$extreme<-sort(out$extreme[!(out$extreme%in%B)])
    }
    out$extreme<-sort(unique(out$extreme[out$extreme>L.base.nodes]))
    out$children<-sort(unique(out$children[out$children>L.base.nodes]))
    out$loglik<-out$modpr<-out$parsRB<-out$pars<-NA
    out
  }
  V.Model.subtract.node<-Vectorize(Model.subtract.node,vectorize.args="node",SIMPLIFY=FALSE)
  
  #Extract the names of models that can be sampled by adding/removing one node to/from model M
  CandidateModels_PG <- function(M){
    M$extreme = (M$extreme)
    M$children = (M$children)
    vec.EnC <- c(M$extreme,M$children)
    vec.MnotE <- M$nodes[!(M$nodes%in%M$extreme)]
    mat.change <- cbind(matrix(0,ncol=length(M$extreme),nrow=length(vec.EnC)),
                        matrix(1,ncol=length(M$children),nrow=length(vec.EnC)))
    diag(mat.change) <- c(rep(1,length(M$extreme)),rep(0,length(M$children)))
    mat.change <- (mat.change+rep(1,length(vec.EnC)))%%2
    modelnames <- apply(mat.change,1,function(r)paste(sort(c(vec.MnotE,vec.EnC[r==1])),collapse=","))
    return(modelnames)
  }
  CandidateModelGlobal <- function(latresp){
    out <- ModelStorage$model[[1]]
    Ch <- out$children
    nodeschange <- Ch[which(rbinom(length(Ch),1,runif(1))==1)]
    while((!is.null(Ch))&(length(nodeschange)>0)){
      for(k in nodeschange){
        out <- Model.add.node(node=k,M=out,SH=SH)
      }
      Ch <- out$children
      nodeschange <- Ch[which(rbinom(length(Ch),1,runif(1))==1)]
    }
    namout=paste(out$nodes,collapse=",")
    nums.instore <- fmatch(namout, names(ModelStorage$model))
    out=ModFit.fn(M=out,latresp)
    
    if(is.na(nums.instore)){
      nums.instore = length(ModelStorage$model)+1
      logpostvec$logpost[[nums.instore]] <- (out$loglik+out$modpr[1])
    }else{
      logpostvec$logpost[[nums.instore]] <- c(unlist(logpostvec$logpost[[nums.instore]]),(out$loglik+out$modpr[1]))
    }
    ModelStorage$model[[namout]] <- out
    out
  }
  #Generate the new models that can be sampled starting from model M
  PG_CandModel <- function(M,SH,latresp){
    new.models <- CandidateModels_PG(M)
    nums.instore <- fmatch(new.models, names(ModelStorage$model))
    not.instore <- new.models[is.na(nums.instore)]
    instore <- new.models[!is.na(nums.instore)]
    
    if(sum(is.na(nums.instore))>0){
      vecall <- c(M$extreme,M$children)
      nodeschange <- vecall[is.na(nums.instore)]
      dropNodes <- M$extreme[M$extreme%in%nodeschange]
      addNodes <- M$children[M$children%in%nodeschange]
      if((length(dropNodes)>0)&(length(addNodes)>0)){
        NewModels=c(V.Model.subtract.node(node=dropNodes,M=M,SH=SH),
                    V.Model.add.node(node=addNodes,M=M,SH=SH))
      }else if(length(addNodes)>0){
        NewModels=V.Model.add.node(node=addNodes,M=M,SH=SH)
      }else if(length(dropNodes)>0){
        dropNodes <- M$extreme[is.na(nums.instore[1:length(M$extreme)])]
        NewModels=V.Model.subtract.node(node=dropNodes,M=M,SH=SH)
      }
      firstpos <- length(ModelStorage$model)+1
      lastpos <- length(ModelStorage$model)+length(NewModels)
      nums.instore[is.na(nums.instore)] <- firstpos:lastpos
      NewModInfo<-VModFit.fn(M=NewModels,latresp)
      ModelStorage$model[firstpos:lastpos] <- NewModInfo
      names(ModelStorage$model)[firstpos:lastpos] = 
        unlist(lapply(NewModels,function(x)paste(x$nodes,collapse=",")))
      
    }
    ModelStorage$model[instore] <- VModFit.fn(M=ModelStorage$model[instore],latresp)
    numpost <- unlist(lapply(ModelStorage$model[nums.instore],function(x)(x$loglik+x$modpr[1])))
    names(numpost)=NULL
    logpostvec$logpost[nums.instore] <- lapply(1:length(nums.instore),function(i){
      c(unlist(logpostvec$logpost[nums.instore[i]]),numpost[i])
    })
    return(list(nums.instore=nums.instore,numerpost=numpost,namesmods=names(ModelStorage$model)[nums.instore] ))
  }
  
  #calculates fwd and bwd proposals and evaluates MH ratio 
  Kernel_PG <- function(M,SH,latresp){ 
    if(runif(1)>shuffle){
      objsfwd <- PG_CandModel(M,SH,latresp)
      postfwd <- exp(objsfwd$numerpost-max(objsfwd$numerpost))
      postfwd <- postfwd/sum(postfwd)
      modsfwd <- objsfwd$nums.instore
      probdrawfwd <- 0.5*(postfwd+(1/length(postfwd)))
      numcand <- sample(modsfwd, 1, replace = FALSE, prob = probdrawfwd)
      Mcand <- ModelStorage$model[[numcand]]
      qfwd <- log(probdrawfwd[modsfwd==numcand])
      objsbwd <- PG_CandModel(Mcand,SH,latresp)
      postbwd <- exp(objsbwd$numerpost-max(objsbwd$numerpost))
      postbwd <- postbwd/sum(postbwd)
      modsbwd <- objsbwd$nums.instore
      probdrawbwd <- 0.5*(postbwd+(1/length(postbwd)))
      qbwd <- log(probdrawbwd[objsbwd$namesmods==paste(M$nodes,collapse=",") ])
      M <- ModelStorage$model[[paste(M$nodes,collapse=",")]]
      margfwd <- Mcand$loglik+Mcand$modpr[1]
      margbwd <- M$loglik+M$modpr[1]
    }else{
      Mcand <- CandidateModelGlobal(latresp)
      margfwd <- Mcand$loglik
      margbwd <- M$loglik
      qbwd = 0
      qfwd = 0
    }
    names(qfwd)=names(qbwd)=names(margfwd)=names(margbwd)=NULL
    if(log(runif(1))< (margfwd - margbwd + qbwd - qfwd) ){
      return(Mcand)
    }else{
      return(M)
    }
  }
}
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#Stochastic Search
{
  model.sequence$modseq[1] <- paste(Model.current$nodes,collapse=",")
  
  function(k,Model.current,lat.response){
      M <- Kernel_PG(M=Model.current,SH=SH,latresp=lat.response)
      Model.current <- M
      model.sequence$modseq[(k+1)] <- paste(Model.current$nodes,collapse=",")
      MatParRB$ParRB[k,] <- Model.current$parsRB
      return(Model.current)
  }
  #return(list(MF.mat=MF.mat,MB.mat=MB.mat,M.start=M.start,nodes=nodes,models=ModelStorage$model,
  #               model.sequence=model.sequence))
}
#-------------------------------------------------------------------------

}
}


