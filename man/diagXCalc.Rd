% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagXCalc.R
\name{diagXCalc}
\alias{diagXCalc}
\title{Calculate the D(Xt) by combining the cell types counts properly.}
\usage{
diagXCalc(xvalues, ncell)
}
\arguments{
\item{xvalues}{Single time-point observation.}

\item{ncell}{The number of cell types over which the process evolves.}
}
\value{
r x r diagonal matrix.
}
\description{
Calculate the D(Xt) by combining the cell types counts properly.
}
\examples{
diagXCalc(xvalues,ncell)
}
