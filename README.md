# Stochastic logistic cell differentiation process - SLCDP

This repository accompanies the manuscript by Pellin et a. __"Tracking hematopoietic stem cell evolution in a Wiskott-Aldrich clinical trial"__.

The abstract of the paper read as follow:

_Hematopoietic Stem Cells (HSCs) are the cells that give rise to all other blood cells and as such, they are crucial in the healthy development of individuals. However, in patients with Wiskott-Aldrich Syndrome (WAS) mutations in the WAS gene lead to a deregulation of hematopoietic cells, resulting in immunodeficiency, increased bleeding, and eczema, as well as increased risks of autoimmune disorders and cancers, such as leukemia and lymphoma. We analyze a revolutionary gene therapy clinical trial, where HSCs harvested from 3 WAS patients' bone marrow have been edited and corrected using viral vectors. Upon re-infusion into the patient, the HSCs multiply and differentiate into other cell types. How this cell multiplication and cell differentiation process operates has until now remained elusive._ 

_This paper aims to model the replenishment of blood lineages resulting from corrected HSCs via a multivariate, density-dependent Markov process and to develop an inferential procedure able to estimate the dynamic parameters given a set of temporally sparsely observed trajectories. Starting from the Master equation, we derive a system of non-linear differential equations for the evolution of the first- and second-order moments over time. We use these moment equations in a generalized method-of-moments framework to perform inference. The performance of the method is evaluated through a simulation study under different sampling scenarios._

_Applying our Markov model and generalized method-of-moments to the Wiskott-Aldrich Syndrome gene therapy clinical trial, resulting in some important results. Contrary to the classical lymphoid/myeloid dichotomy theory of stem cell differentiation, our method supports a recent myeloid-based model that claims that the fates of lymphoid and myeloid cells remain coupled after the loss of erythroid potential._

## How to install

The package SLCDP requires the installation of [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio), software available at no-cost under the [IBM academic initiative](https://www.ibm.com/academic/home).
CPLEX is called by R through the dedicated `R` package `Rcplex`.
For a smooth installation, we suggest following the following steps:
1. Register to IBM academic initiative.
2. Go to the [search page](https://www-03.ibm.com/isc/esd/dswdown/home.wss) and search the following _part number_: CRW5CML 
3. Download the IBM ILOG CPLEX Optimization Studio V12.6.3 Multiplatform Multilingual eAssembly (CRW5CML) that fits your system. Newer versions (>12.x) are not properly recognized by Rcplex.
4. Install IBM ILOG CPLEX Optimization Studio following the instruction. Take note of the installation directory.
5. Download `Rcplex` package from [CRAN](https://cran.r-project.org/web/packages/Rcplex/index.html) and install it using the command line:
`R CMD INSTALL --configure-args="--with-cplex-dir=CPLEX_INSTALLATION_DIRECTORY" Rcplex_X.X.tar.gz`
(CPLEX_INSTALLATION_DIRECTORY is /opt/ibm/ILOG/CPLEX_Studio_Preview126/cplex using default setting on linux machine)
6. Test the installation by running the example included in `Rcplex` help function after having started a new `R` session and loaded `Rcplex` library.
7. Install `SLCDP` package using:
```
library(devtools)
install_github("dp3ll1n/SLCDP")
```
8. Run the example available in the Vignette
