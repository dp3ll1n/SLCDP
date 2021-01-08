#' Calculate the solution of the ODEs
#'
#' Solve the system of differential equations given initial values and parameters.
#' @param initialValues1_2ODEs Differential equations initial values.
#' @param a ODEs parameters values.

#' @return  List of prediction.
#' @export
#' @examples
#' calculatePred(initialValues1_2ODEs,a)


calculatePred=function(initialValues1_2ODEs,a){
  sdePred=lapply(split(initialValues1_2ODEs,seq(nrow(initialValues1_2ODEs))),function(x) {
    (as.numeric(odeintr::integrate_sys(dxdt,as.numeric(x),dt, step_size = stepSize, start = 0, observer = named_rec,atol=toll,rtol=toll)[-1]))})
  return(sdePred)
}
