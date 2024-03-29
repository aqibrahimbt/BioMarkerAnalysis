% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_analysis.R
\name{dds_heatmap_gene_counts}
\alias{dds_heatmap_gene_counts}
\title{Generates heatmap for top {n} expressed genes and outputs to JSON file
based on normalized count data}
\usage{
dds_heatmap_gene_counts(dds, subset, datasets, outdir)
}
\arguments{
\item{dds}{- dds object}

\item{subset}{- number of genes required}

\item{datasets}{- datasets used for the analysis}

\item{outdir}{- output directory}
}
\description{
Generates heatmap for top {n} expressed genes and outputs to JSON file
based on normalized count data
}
