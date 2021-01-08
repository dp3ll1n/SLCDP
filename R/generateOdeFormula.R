#' Generate SLCDP ODE formulas for moments evolution. 
#'
#' Automatically generate the ODE formulas describing the first and second-order moments evolution for a SLCDP with ncell cell types. 
#' @param ncell The number of cell types over which the process evolves.
#' @return Text corresponding to formula.
#' @export
#' @examples
#' generateOdeFormula(5)


generateOdeFormula=function(ncell){
  allr=rep(0,ncell)
  for (i in c(1:ncell)){
    for (a in c(1:ncell)) {       
      r=rep(0,ncell)
      if(i!=a){r[i]=1;r[a]=-1}       
      allr=cbind(allr,r)
    }
  }
  allr=allr[,-1]
  rm(r)
  
  duM=diag(1,ncell,ncell)
  deM=diag(-1,ncell,ncell)
  netMat=t(cbind(duM,deM,allr))
  
  param=paste("a",1:((ncell*(ncell+2))),sep="")
  deriv=paste("dX",1:ncell,sep="")
  proc=paste("X",1:ncell,sep="")
  xrep=rep(proc,(ncell+2))
  
  eqm=c()
  eq=""
  
  for (i in 1:dim(netMat)[2]){
    eq=""
    for (a in 1:dim(netMat)[1]){
      if(a %in% c((ncell+1):(ncell*2))){if(netMat[a,i]!=0) {eq=paste(eq,paste("(",paste(netMat[a,i],param[a],paste(" ",xrep[a],xrep[a]," ",sep=""),sep="*"),")",sep=""),sep="+")} }
      else {if(netMat[a,i]!=0) { eq=paste(eq,paste("(",paste(netMat[a,i],param[a],paste(" ",xrep[a]," ",sep=""),sep="*"),")",sep=""),sep="+")}}
    }
    eqm[i]=eq   
  }
  eqall=eqm
  ###
  eqm=matrix(0,ncell,ncell)
  eq=""
  for (j in 1:length(proc)){
    for (i in 1:dim(netMat)[2]){
      eq=""
      for (a in 1:dim(netMat)[1]){
        if(a %in% c((ncell+1):(ncell*2))){if(netMat[a,i]!=0) {eq=paste(eq,paste("(",paste(netMat[a,i],param[a],paste("( ",proc[j]," *( ",xrep[a],"",xrep[a]," - ",xrep[a]," ^2)+2* ",xrep[a]," *( ",xrep[a],"",proc[j]," -( ",xrep[a]," * ",proc[j]," ))+(( ",xrep[a]," ^2)* ",proc[j]," ))",sep=""),sep="*"),")",sep=""),sep="+")}
        }
        else {if(netMat[a,i]!=0) {eq=paste(eq,paste("(",paste(netMat[a,i],param[a],paste(" ",proc[j],xrep[a]," ",sep=""),sep="*"),")",sep=""),sep="+")}}
      }
      eqm[j,i]=paste(eq,sep="=")   
    }
  }
  
  ball=c()
  
  for (i in c(1:ncell)){
    for (a in c(1:ncell)){
      b1=eqm[i,a]
      b2=eqm[a,i]
      b=paste(b1,b2,sep="")   
      ball=c(ball,b)
    }
    
  }
  
  eqm=matrix(0,ncell,ncell)
  eq=""
  for (j in 1:dim(netMat)[2]){
    for (i in 1:dim(netMat)[2]){
      eq=""
      for (a in 1:dim(netMat)[1]){
        if(a %in% c((ncell+1):(ncell*2))){if(netMat[a,i]!=0 && netMat[a,j]!=0) { eq=paste(eq,paste("(",paste(netMat[a,i],param[a],paste(" ",xrep[a],xrep[a]," ",sep=""),netMat[a,j],sep="*"),")",sep=""),sep="+")} }
        else {if(netMat[a,i]!=0 && netMat[a,j]!=0) { eq=paste(eq,paste("(",paste(netMat[a,i],param[a],paste(" ",xrep[a]," ",sep=""),netMat[a,j],sep="*"),")",sep=""),sep="+") } }
      }
      eqm[j,i]=eq
    }
  }
  
  eq2momentProv=paste(ball,as.vector(eqm),sep="")
  
  for (i in c(1:ncell)){
    for (a in c((i+1):ncell)){
      eq2momentProv=gsub(as.character(paste(" ",proc[a],proc[i]," ",sep="")),as.character(paste(" ",proc[i],proc[a]," ",sep="")), eq2momentProv)
    }
    
  }
  
  eqmVarCov=matrix(eq2momentProv,ncell,ncell,byrow=T)
  eqmVarCovUT=as.vector(eqmVarCov[lower.tri(eqmVarCov,diag=T)])
  eqall=c(eqall,eqmVarCovUT)
  
  #############
  
  comb=expand.grid(1:ncell,1:ncell)
  ind=apply(comb,1,function(x) ifelse(x[2]<=x[1],1,0))
  combOk=comb[which(ind==1),]
  eqName=paste("X",1:ncell,sep="")
  eqName=c(eqName,apply(combOk,1,function(x) {paste0(eqName[x[2]],eqName[x[1]])}))
  
  
  elem=as.vector(eqName)
  subelem=paste("x[",1:(length(elem)),"]",sep="")
  param=paste("a",1:((ncell*(ncell+2))),sep="")
  subparam=paste("a[",1:(length(param)),"]",sep="")
  
  expelem=elem[1:ncell]
  expelem=paste(" ",expelem," ","\\^2",sep="")
  subexpelem=paste(" ",paste(elem[1:ncell]," * ",sep=""),elem[1:ncell]," ",sep="") 
  
  form=data.frame(V1=eqall)
  
  for(i in 1:dim(form)[1]){
    ooo=as.character(form[i,1])
    for(a in length(expelem):1){
      ppp=gsub(expelem[a],subexpelem[a],(ooo))
      ooo=ppp
    }
    form[i,1]=as.character(ooo)
  }
  
  
  for(i in 1:dim(form)[1]){
    ooo=as.character(form[i,1])
    for(a in length(elem):1){
      ppp=gsub(elem[a],subelem[a],(ooo))
      ooo=ppp
    }
    form[i,1]=as.character(ooo)
  }
  
  for(i in 1:dim(form)[1]){
    ooo=as.character(form[i,1])
    for(a in length(param):1){
      ppp=gsub(param[a],subparam[a],(ooo))
      ooo=ppp
    }
    form[i,1]=as.character(ooo)
  }
  
  form=apply(form,1,function(x) gsub(" ", "", x, fixed = TRUE))
  formPaste=paste(form,collapse=",")
  finform=paste("dxdt = function(x, t) c(", formPaste, ")",sep="",collapse = "")
  hh=parse(text=finform)
  return(hh)
}
