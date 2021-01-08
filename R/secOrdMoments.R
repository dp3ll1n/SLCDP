#' Combination for second order moments 
#'
#' @param ncell The number of cell types over which the process evolves.

#' @return  Matrix.
#' @export
#' @examples
#' secOrdMoments(ncell)


secOrdMoments=function(ncell){
  comb=expand.grid(1:ncell,1:ncell)
  ind=apply(comb,1,function(x) ifelse(x[2]<=x[1],1,0))
  combOk=comb[which(ind==1),]
  return(combOk)
}
