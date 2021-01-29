#' Create a list containing rates constraints and bounds for the definition of the Quadratic Programming problem in the Cplex optimization environment.
#'
#' It initializes and stores constraints and bounds for each individual rates. By default, given ncell, it sets all reates to be non-negative, exept for differentiation rates matrix main diagonal elements set to 0. Using this function it is possible to set specific rates to fixed values (non-negative).
#' @param ncell The number of cell types over which the process evolves.
#' @param duplRatesLB duplRatesUB deathRatesLB deathRatesUB Vectors of length ncell containing duplication (dupl) and death (death) rates upper- (UB) and lower- (LB) bounds.
#' @param diffRatesLB diffRatesUB Matrix of size ncell x ncell containing differentiation rates upper- (UB) and lower- (LB) bounds.
#' @return  List with A, (inequality constraints),lb (lower bound), ub (upper bound).

#' @export
#' @examples
#' setConstraintsOptiSolve(ncell,duplRatesLB=NA,duplRatesUB=NA,deathRatesLB=NA,deathRatesUB=NA,diffRatesLB=NA,diffRatesUB=NA)


setConstraintsOptiSolve=function(ncell,
                        duplRatesLB=NA,
                        duplRatesUB=NA,
                        deathRatesLB=NA,
                        deathRatesUB=NA,
                        diffRatesLB=NA,
                        diffRatesUB=NA){
  if(is.na(duplRatesLB)){
    duplRatesLB=rep(0,ncell)
  }
  if(is.na(duplRatesUB)){
    duplRatesUB=rep(NA,ncell)
  }
  if(is.na(deathRatesLB)){
    deathRatesLB=rep(0,ncell)
  }
  if(is.na(deathRatesUB)){
    deathRatesUB=rep(NA,ncell)
  }
  if(is.na(diffRatesLB)){
    diffRatesLB=matrix(0,nrow = ncell,ncol = ncell)
  }
  if(is.na(diffRatesUB)){
    diffRatesUB=matrix(NA,nrow = ncell,ncol = ncell)
  }
  if(!(all(c(length(duplRatesLB),
             length(duplRatesUB),
             length(deathRatesLB),
             length(deathRatesUB))==ncell))){
    stop("Duplication and death rates upper and lower bounds must have ncell length")
  }

  if(!(all(c(ncol(diffRatesLB),
             nrow(diffRatesLB),
             ncol(diffRatesUB),
             nrow(diffRatesUB))==ncell))){
    stop("Differentiation rates upper and lower bounds must be a matrices of size ncell x ncell")
  }
  ### constrain matrix QP ###
  diag(diffRatesUB)=0
  A=diag(ncell*(ncell+2))
  lb=c(duplRatesLB,deathRatesLB,as.vector(diffRatesLB))
  ub=c(duplRatesUB,deathRatesUB,as.vector(diffRatesUB))
  constraints=list(A=A,lb=lb,ub=ub)
  return(constraints)
}
