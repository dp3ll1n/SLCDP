#' Extract clones observations at time points t= 1 .. n. 
#'
#' Extract clones observations at time points t= 1 .. n. These are the values to be matched by ODEs solutions. 
#' @param clonesTrajs Set of clone trajectories.
#' @return  Matrix.
#' @export
#' @examples
#' endValuesColVec(clonesTrajs)


endValuesColVec=function(clonesTrajs){
  clonesTrajsEndValues=do.call(rbind,lapply(clonesTrajs,function(x) {x[-1,]} ))
  clonesTrajsEndValues=as.vector(t(as.matrix(clonesTrajsEndValues)))
  return(clonesTrajsEndValues)
}
