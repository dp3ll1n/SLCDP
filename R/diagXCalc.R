#' Calculate the D(Xt) by combining the cell types counts properly. 
#'
#' @param xvalues Single time-point observation.
#' @param ncell The number of cell types over which the process evolves.

#' @return  r x r diagonal matrix.
#' @export
#' @examples
#' diagXCalc(xvalues,ncell)


diagXCalc=function(xvalues,ncell){
  diagX=diag( c(xvalues,xvalues^2,rep(xvalues,(ncell))) )
  return(diagX)
}
