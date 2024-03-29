% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workbench_advise_lipidomics.R
\name{workbench_advise_lipidomics}
\alias{workbench_advise_lipidomics}
\title{workbench_advise_lipidomics}
\usage{
workbench_advise_lipidomics(filtered, coldata, metadata, data_type, id)
}
\arguments{
\item{filtered}{List. It is the result from the \code{na_advise_lipidomics} function.}

\item{coldata}{the colData object from the Metabolomics Workbench output.}

\item{metadata}{metadata to be saved in the final summarized experiment.}

\item{data_type}{information about the data type (Area, Concentration, Peak Heigth...).}

\item{id}{the id of the selected metabolomics workbench experiment .}
}
\value{
res: a list with two SummerizedExperiment
}
\description{
\code{workbench_advise_lipidomics} is the function for the creation of the
SummarizedExperiment object from a Metabolomics Workbench output.
}
\details{
This function take the filtered matrix (assay), the target file (coldata) and 
the metadata to build a SummarizedExperiment. The rowData of the SumExp is also builded.
}
\note{
Last change 28/12/2021
}
