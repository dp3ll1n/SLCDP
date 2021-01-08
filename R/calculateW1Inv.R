#' Calculate the global block diagonale covariance matrix inverse from ODEs solutions
#'
#' Combines differential equation solutions and calculate the inverse of covariance matrix W. It is possible to add additional variance on the main diagonal (addtional source of noise) by setting lambda.
#' @param sdePred Differential equations solutions.
#' @param lambda Additional variance to be added on the main diagonal of the global covariance matrix.
#' @param ncell The number of cell types over which the process evolves.
#' @param combOk The number of cell types over which the process evolves.

#' @return  Block diagonal matrix.
#' @export
#' @examples
#' calculateW1Inv(sdePred,lambda=0,ncell,combOk)


calculateW1Inv=function(sdePred,lambda=0,ncell,combOk){
  w=lapply(sdePred,function(x) {calcW(x,ncell,combOk)})
  w1Inv=lapply(w,function(x){diag(x)=diag(x)+lambda; y=(ginv(x)); return(y)} )
  w1Inv=(Matrix::bdiag(w1Inv))
  return(w1Inv)
}
