#' Calculate the starting values for moments ODEs.
#'
#'
#' @param clonesTrajsStartValues Clone trajectories starting value returned by startValues()
#' @param ncell The number of cell types over which the process evolves.

#' @export
#' @examples
#' initialValues1_2ODEsCalc(clonesTrajsStartValues,ncell)


initialValues1_2ODEsCalc=function(clonesTrajsStartValues,ncell){
  comb=expand.grid(1:ncell,1:ncell)
  ind=apply(comb,1,function(x) ifelse(x[2]<=x[1],1,0))
  combOk=comb[which(ind==1),]
  eqName=paste("X",1:ncell,sep="")
  eqName=c(eqName,apply(combOk,1,function(x) {paste0(eqName[x[2]],eqName[x[1]])}))
  initStateDF=lapply(split(clonesTrajsStartValues,seq(nrow(clonesTrajsStartValues))),function(y) { as.numeric(unlist(apply(combOk,1,function(x,y) {y[x[1]]*y[x[2]]},y)))})
  initStateDF=do.call(rbind,initStateDF)
  initStateDF=cbind(clonesTrajsStartValues,initStateDF)
  colnames(initStateDF)=eqName
  return(initStateDF)
}
