% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.r
\name{get.data}
\alias{get.data}
\title{Get transcript count matrix, conditions and feature length of requested dataset}
\usage{
get.data(datasets, meta_data, working_dir, file = "counts.1.tsv",
  formula = as.formula("~condition"), transcript_selection = NA,
  filter = NA, filter_column = "disease", pseudo_count = 0,
  min_counts = 1)
}
\arguments{
\item{meta_data}{meta data table}

\item{working_dir}{base folder including all datasets}

\item{file}{name of transcript count file}

\item{formula}{formula to select conditions}

\item{transcript_selection}{transcript selection pattern to grep specific transcripts}

\item{filter}{filter data}

\item{filter_column}{column to grep by filter}

\item{pseudo_count}{add pseudo count to feature counts}

\item{min_counts}{minimum sum of counts per feature}

\item{dataset}{dataset numbers}
}
\description{
Get transcript count matrix, conditions and feature length of requested dataset
}
