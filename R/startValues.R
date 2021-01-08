#' Extract clones observations at time points t= 0 .. (n-1). 
#'
#' Extract clones observations at time points t= 0 .. (n-1), useful for the calculation of ODEs starting values. 
#' @param clonesTrajs Set of clone trajectories.
#' @return  Matrix.
#' @export
#' @examples
#' startValues(clonesTrajs)


startValues=function(clonesTrajs){
  clonesTrajsStartValues=do.call(rbind,lapply(clonesTrajs,function(x) {x[-dim(x)[1],]} ))
  return(clonesTrajsStartValues)
}
