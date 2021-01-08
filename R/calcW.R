#' Calculate covariance matrix
#'
#' Combines differential equation solutions as encoded in combOk matrix and calculate the covariance matrix W
#' @param diffEqSol Differential equations solutions.#'
#' @param ncell The number of cell types over which the process evolves.
#' @param combOk The number of cell types over which the process evolves.

#' @return  Matrix.
#' @export
#' @examples
#' calcW(diffEqSol,ncell,combOk)


calcW=function(diffEqSol,ncell,combOk){
  moment1=diffEqSol[1:ncell]
  moment2=diffEqSol[(ncell+1):length(diffEqSol)]
  m1_m1=apply(combOk,1,function(x,moment1) {moment1[x[1]]* moment1[x[2]]},moment1)
  wi=diag(ncell)
  wi[lower.tri(wi, diag=TRUE)] <- (moment2-m1_m1)
  wi[upper.tri(wi)] <- t(wi)[upper.tri(wi)]
  return(wi)	
}
