% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{rce_step}
\alias{rce_step}
\title{RCE step to remove clusters with the lowest accuracy.}
\usage{
rce_step(svm_result, kclusters, train)
}
\arguments{
\item{svm_result}{result from svm scoring step}

\item{kclusters}{required number of clusters}

\item{train}{training dataset}
}
\description{
RCE step to remove clusters with the lowest accuracy.
}
