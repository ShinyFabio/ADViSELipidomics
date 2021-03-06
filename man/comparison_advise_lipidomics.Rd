% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparison_advise_lipidomics.R
\name{comparison_advise_lipidomics}
\alias{comparison_advise_lipidomics}
\title{comparison_advise_lipidomics}
\usage{
comparison_advise_lipidomics(out, gradient = FALSE)
}
\arguments{
\item{out}{List. It is the result from the \code{diffexp_advise_lipidomics} function.}

\item{gradient}{Logical value. If TRUE, the Venn Diagram is colored taking into account
the number of common lipids by a gradient. Default = FALSE.}
}
\value{
Venn diagram and Upset plot of common lipids among different comparisons
(comparison of all the selected contrasts).
}
\description{
\code{comparison_advise_lipidomics}
}
\details{
The Venn Diagramm can be in classical form or in gradient form.
}
\note{
Last change 17/12/2021
}
