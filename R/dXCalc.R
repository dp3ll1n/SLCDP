#' Calculate the differences in cell counts among consecutive timepoints. 
#'
#' Calculate the differences in cell counts among consecutive observations for a set of clone trajectories. The output is in the OLS procedure to set parameters initial values. 
#' @param clonesTrajs Set of clone trajectories.
#' @return  Deltas matrix.
#' @export
#' @examples
#' dXCalc(clonesTrajs)


dXCalc=function(clonesTrajs){
  dX=lapply(clonesTrajs,function(x) diff(as.matrix(x)))
  dX=unlist(lapply(dX,function(x) as.vector(t(x))))
  return(dX)
}
