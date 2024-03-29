% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.r
\name{find.collinearity}
\alias{find.collinearity}
\title{Builds a model matrix of the condition matrix and the formula. 
Tests if the Variance inflation factor of each factor is below threshold. Otherwise, 
the factor is collinear to other fator. Remove factor from condition matrix.}
\usage{
find.collinearity(coldata, formula, maxVif = 10)
}
\arguments{
\item{coldata}{condition matrix of dataset}

\item{formula}{formula}

\item{maxVif}{VIF threshold}
}
\description{
Builds a model matrix of the condition matrix and the formula. 
Tests if the Variance inflation factor of each factor is below threshold. Otherwise, 
the factor is collinear to other fator. Remove factor from condition matrix.
}
