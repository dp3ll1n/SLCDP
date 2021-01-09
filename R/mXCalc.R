#' Calculate the linear predictor matrix M for the local linear approximation.
#'
#' By combining cell count observations, net effect matrix and dt it returns the matrix of linear predictors for setting parameters initial values.
#' @param clonesTrajs Set of clone trajectories.
#' @param netMat Net effect matrix.
#' @param dt Time between observation
#' @param ncell The number of cell types over which the process evolves.

#' @return  Predictors matrix.
#' @export
#' @examples
#' mXCalc(clonesTrajs,netMat,dt,ncell)


mXCalc=function(clonesTrajs,netMat,dt,ncell){
  MXols=lapply(split(clonesTrajs,seq(nrow(clonesTrajs))),function(x,ncell) {(t(netMat)) %*% diagXCalc(x,ncell) *dt},ncell)
  MXols <- do.call(rbind,MXols)
  MXols=Matrix::Matrix(MXols,sparse=T)
  return(as.matrix(MXols))
}
