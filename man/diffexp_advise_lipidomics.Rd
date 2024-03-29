% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffexp_advise_lipidomics.R
\name{diffexp_advise_lipidomics}
\alias{diffexp_advise_lipidomics}
\title{diffexp_advise_lipidomics}
\usage{
diffexp_advise_lipidomics(
  out,
  rep_mean = FALSE,
  rep_effect = TRUE,
  bs_norm = "none",
  batch_type = "remove",
  batch_method = "limma",
  batch_vars = "",
  eb_robust = FALSE,
  eb_trend = FALSE,
  thresh = 0.05,
  decide_met = "separate"
)
}
\arguments{
\item{out}{List. It is the result from the \code{expdesign_advise_lipidomics} function.}

\item{rep_mean}{Logical value. It is for the selection of the dataset with the replicates
(FALSE) or the dataset with the samples (TRUE). Default = FALSE.}

\item{rep_effect}{Logical value. If rep_mean == FALSE, it is possible to consider the
influence of the replicates presence for fitting the linear model, taking into account a
value of correlation among the replicates. Default = TRUE.}

\item{bs_norm}{Character. It is possible to normalize data between replicates
or samples. Possible options are "none", "scale" and "quantile" (only if it isn't a concentration matrix). 
Default = "none".}

\item{batch_type}{Character string. Different methodologies to cope with the eventual presence
of batch effects: "remove", "fit". Default = "remove".}

\item{batch_method}{Algorithm used for the batch effect. Choices are: "limma", "combat_nonparam" and "combat_param".}

\item{batch_vars}{Character vector, with the names of the variables from the target file,
considered as batch effects for the entire experiment. Default = "".}

\item{eb_robust}{Logical value for inner use of empirical Bayes moderation of the
standard errors towards a global value. More information in "limma' package
(TRUE for RNAseq data). Default = FALSE.}

\item{eb_trend}{Logical value for inner use of empirical Bayes moderation of the
standard errors towards a global value. More information in "limma' package
(TRUE for RNAseq data). Default = FALSE.}

\item{thresh}{Numerical value. Threshold on adjust P-value in order to consider a lipid
differentially expressed or not. Default = 0.05.}

\item{decide_met}{Character string specifying how lipids and contrasts are to be combined
in the multiple testing scheme. Choices are "separate", "global", "hierarchical" or "nestedF".
More information in "limma' package. Default = "separate".}
}
\value{
res: a list with results from the experimental design, updated with the results
from fitting the linear model and appling decisional test on the results from the model.
Moreover, all the tables and plots for differential expressed lipids are available.

res: a list with results from the experimental design, updated with the results
from fitting the linear model and appling decisional test on the results from the model.
Moreover, all the tables and plots for differential expressed lipids are available.
}
\description{
\code{diffexp_advise_lipidomics} is the function for creation of the experimental
design.

\code{diffexp_advise_lipidomics} is the function for eventually removing the
batch effects and applying the differential analysis per lipids among replicates/samples,
taking into account the comparisons of interest.
}
\details{
This function has many tasks. At first, it log-transforms data for subsequent analysis
and provides a simple method for normalizing data (scaling) per replicates or samples.
Then, it is possible to remove batch effects by "limma" method (at most two batch variables) or
by ComBat method (from "sva" package), in parametric or non-parametric modes. In particular
the option "combat_param_full" will run ComBat with parametric adjustment and the null model
with design matrix, "combat_nonparam_full" will run ComBat with nonparametric adjustment and
the null model with design matrix (non-parametric approach can be very time consuming).
After this correction, the linear model is fitted and the result are reported in form of table
(raw tables and filtered tables) and in form of plots (MA-plots, Volcano plots, number of
differential expressed lipids plots), separated for all the comparisons of interest.
}
\note{
last change 17/12/2021
}
