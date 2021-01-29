## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE---------------------------------------------------------
write_matex <- function(x) {
  begin <- "$$\\begin{bmatrix}"
  end <- "\\end{bmatrix}$$"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  writeLines(c(begin, X, end))
}

## ----setup--------------------------------------------------------------------
library(SLCDP)
ncell=5

diffRate=matrix( c(0,0.2,0.35,0,0,0,0,0,0.75,0,0,0,0,0.25,0.5,0,0,0,0,0,0,0,0,0,0),ncell,ncell,byrow = T)
dupRate =c(1.000,1.500,1.800,2.500,2.800)
deathRate=c(0.033,0.030,0.045,0.031,0.043);
E_0=rep(0,ncell)
E_0[1]=1

tEnd=10	# time end
dt=1	# delta time

## ----results='asis'-----------------------------------------------------------
eps=1e-4
toll=1e-3
stepSize=dt/10
CombOk=secOrdMoments(ncell)
form=generateOdeFormula(ncell)
eval(form)

## ----results='asis'-----------------------------------------------------------
NetEffectMatrix= netEffectMatrixCalc(ncell)

## -----------------------------------------------------------------------------
generateClone(ncell,NetEffectMatrix, dupRate, deathRate, diffRate, dt, tEnd, E_0)

## ----results='asis'-----------------------------------------------------------
N=100
ClonesTrajs=lapply(c(1:N), function(x) { 
  generateClone(ncell,NetEffectMatrix, dupRate, deathRate, diffRate, dt, tEnd, E_0)
  })

## ----results='asis'-----------------------------------------------------------
ConstraintsCplex=setConstraintsCplex(5)

## ----results='asis'-----------------------------------------------------------
deltaX=dXCalc(ClonesTrajs)
ClonesTrajsStartValues=startValues(ClonesTrajs)
ClonesTrajsEndValues=endValuesColVec(ClonesTrajs)
MXCalc=mXCalc(ClonesTrajsStartValues,NetEffectMatrix,dt,ncell)
InitialValues1_2ODEs=initialValues1_2ODEsCalc(ClonesTrajsStartValues,ncell)

## ----results='asis'-----------------------------------------------------------
a=runOptimizationCplex(MXCalc,deltaX,W1Inv=NULL,ConstraintsCplex)

## ----results='asis'-----------------------------------------------------------
SdePred=calculatePred(InitialValues1_2ODEs,a)
W1=calculateW1Inv(SdePred,lambda=0,ncell,CombOk)
#a=runOptimization(MXCalc,ClonesTrajsEndValues,W1Inv=W1,Constraints)
a=runOptimizationCplex(MXCalc,deltaX,W1Inv=W1,ConstraintsCplex)

## ----results='asis'-----------------------------------------------------------
SdePred=calculatePred(InitialValues1_2ODEs,a)
W1=calculateW1Inv(SdePred,lambda=0,ncell,CombOk)
FX=calculateCountsPred(SdePred,ncell)
aCurr=a
iter=0
maxiter=100
while( (sum(abs(a-aCurr))>0.05 && (iter<maxiter))|( iter==0) ){
  iter=iter+1
  print(sum(abs(a-aCurr)))
  aCurr=a
  J=matrix(0,nrow = length(FX),ncol = length(aCurr))
  for(parj in which(aCurr!=0)){
    a=aCurr
    a[parj]=aCurr[parj]+eps
    JtP=calculatePred(InitialValues1_2ODEs,a)
    JtPv=calculateCountsPred(JtP,ncell)
    J[,parj]=(JtPv-FX)/(eps)
  }
  ConstraintsCplex$lb=-a
  da=runOptimizationCplex(J,(ClonesTrajsEndValues-FX),W1Inv=W1,ConstraintsCplex,bvec=a)
  a=a+1*(da)
  SdePred=calculatePred(InitialValues1_2ODEs,a)
  W1=calculateW1Inv(SdePred,lambda=0,ncell,CombOk)
  FX=calculateCountsPred(SdePred,ncell)
  print(a)

}


