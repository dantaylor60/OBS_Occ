require(combinat)
allterm.generator <-function(v,r){
  do.call(rbind,lapply(1:r,FUN=function(v,x){
    t(xsimplex(p=v,n=x))},v=v))
}

nodes <- function(x,O1.names){
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

node.edges <- function(x,polymat,vertices,degree.vec,degree.max){
  if(degree.vec[x]< degree.max){
    do.call(rbind,lapply(1:dim(polymat)[2],FUN=function(u){
      pos.x = which(vertices%in%x)
      poly.xto = polymat[pos.x,]
      poly.xto[u] = poly.xto[u]+1
      c(x,names(which(apply(polymat,MARGIN=1,identical,y=poly.xto))))
    }))
  }
}

MF.DAG <- function(X.linear,maxdeg,exclusions=NA,MF.formula=NA){
  deg1names = names(X.linear)
  var.class = c(lapply(1:length(deg1names),function(x)class(X.linear[,x])),recursive=T)
  categorical = deg1names[var.class=="factor"]
  if(!identical(categorical,character(0)) && maxdeg>1){ 
    categorical = paste(categorical,"^2",sep="")
    exclusions = c(categorical,exclusions)
    exclusions = unique(exclusions[!is.na(exclusions)])
  }
  nvars = length(deg1names)
  polyterms = rbind(rep(0,nvars),allterm.generator(v=nvars,r=maxdeg))
  vertices = apply(polyterms,MARGIN=1,FUN = nodes, O1.names=deg1names )
  form.names = vertices["form.alias",]
  vertices = vertices["nam",]
  rownames(polyterms) = vertices 
  degree = rowSums(polyterms)
  edges = do.call(rbind,lapply(vertices,FUN = node.edges,polymat=polyterms,
                               vertices=vertices,degree.vec=degree,degree.max=maxdeg))
  if(sum(is.na(exclusions))==0){
    
    mindeg.excl = min(degree[vertices%in%exclusions])
    
    for(i in mindeg.excl:(maxdeg+1)){
      rmv = do.call(rbind,lapply(exclusions,FUN=node.edges,
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


model.enumerate<-function(MB,MF){
  
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

