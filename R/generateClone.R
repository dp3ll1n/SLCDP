#' @useDynLib SLCDP
#' @export
#'
#'

generateClone=function(NCell, M, dupRate, deathRate, diffRate, dt, tEnd, E0){
  generateCloneTraj(NCell, M, dupRate, deathRate, diffRate, dt, tEnd, E0)
}
