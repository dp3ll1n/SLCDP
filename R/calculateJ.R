#' Calculate the Jacobian matrix
#'
#' @param initialValues1_2ODEs Differential equations initial values.
#' @param a ODEs parameters values.
#' @param ncell The number of cell types over which the process evolves.
#' @param eps Delta for numeric derivatives calculation 

#' @return  Jacobian matrix.
#' @export
#' @examples
#' calculateJ(initialValues1_2ODEs,a,ncell,eps)


calculateJ=function(initialValues1_2ODEs,a,ncell,eps){
  ainit=a
  J=NULL
  for(parj in which(ainit!=0)[1]){
    a=ainit
    a[parj]=ainit[parj]+eps
    JtP=calculatePred(initialValues1_2ODEs,a)
    JtPv=calculateCountsPred(JtP,ncell)
    a[parj]=ainit[parj]-eps
    JtM=calculatePred(initialValues1_2ODEs,a)
    JtMv=calculateCountsPred(JtM,ncell)  
    J=(JtPv-JtMv)/(2*eps)  
  }
  a=ainit
  return(JtMv)
}
