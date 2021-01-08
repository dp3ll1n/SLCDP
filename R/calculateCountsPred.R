#' Extract first order moments values from ODEs solutions 
#'
#' @param sdePred Differential equations solutions.
#' @param ncell The number of cell types over which the process evolves.

#' @return  Vector of predicted cell counts.
#' @export
#' @examples
#' calculateCountsPred(sdePred,ncell)


calculateCountsPred=function(sdePred,ncell){
  fX=lapply(sdePred, function(x) x[1:ncell])
  fX=unlist(lapply(fX,function(x) as.vector(t(x))))
  return(fX)
}
