#' Generate the net effect matrix for the SLCDP. 
#'
#' Generate the r x ncell net effect matrix describing the change in cell counts induces by each of the r possible cellular events. The first ncell rows correspond to duplication events (+1 in process state), followed by ncell death events (-1) and all possible differentiation (-1;+1 respectively for the cell type of origin and destination). (). , the  ODE formulas describing the first and second-order moments evolution for a SLCDP with ncell cell types. 
#' @param ncell The number of cell types over which the process evolves.
#' @return  r x ncell matrix.
#' @export
#' @examples
#' netEffectMatrixCalc(5)


netEffectMatrixCalc=function(ncell){
  allr=rep(0,ncell)
  for (i in c(1:ncell)){
    for (a in c(1:ncell)) {		
      r=rep(0,ncell)
      if(i!=a){r[i]=1;r[a]=-1}		
      allr=cbind(allr,r)
    }
  } 
  allr=allr[,-1]
  rm(r)
  
  duM=diag(1,ncell,ncell)
  deM=diag(-1,ncell,ncell)
  netMat=t(cbind(duM,deM,allr))
  return(netMat)
}
