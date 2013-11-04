\name{strumFittedModel-class}
\docType{class}

\alias{strumFittedModel-class}
\alias{fittedParametersCovMatrix<-}
\alias{fittedParametersCovMatrix}
\alias{fittedParameters<-}
\alias{fittedParameters}
\alias{fittedParametersCovMatrix<-,strumFittedModel-method}
\alias{fittedParametersCovMatrix,strumFittedModel-method}
\alias{fittedParameters<-,strumFittedModel-method}
\alias{fittedParameters,strumFittedModel-method}
\alias{show,strumFittedModel-method}

\title{Class \code{"strumFittedModel"}}
\description{
\code{strumSimModel} is an S4 class that contains the result from a strum analysis, i.e., information of a trait model with fitted parameter values.
}
\section{Objects from the Class}{
Objects should be created by the calls to the \code{\link{strum}} function.
}
\section{Slots}{
  \describe{
    \item{\code{myStrumModel}:}{Object of class \code{"strumModel"} stating the analysis model.}
    \item{\code{modelValidity}:}{Indicator stating the validity of the analysis model.}
    \item{\code{fittedParameters}:}{Object of class \code{"data.frame"} containing the fitted parameter values of the model with the standard errors, confidence intervals, and p-values.}
    \item{\code{fittedParametersCovMatrix}:}{The variance-covariance matrix of the fitted parameters.}
    \item{\code{deltaParameters}:}{For internal use}
    \item{\code{deltaParametersCovMatrix}:}{For internal use}
    \item{\code{parDiff}:}{For internal use}
    \item{\code{parDiffCovMatrix}:}{For internal use}
    \item{\code{chiTestOut}:}{Object of class \code{"data.frame"} containing the un-adjusted(1), mean scaled(2), and mean/variance scaled(3) chi-square index of fits with the degrees of freedom and p-values along with the exact p-value(4) from simulation.}
  }
}
\section{Methods}{
  \describe{
    \item{fittedParameters}{\code{signature(object = "strumFittedModel")}: Accessor function, returns the fittedParameters of the object.}
    \item{fittedParametersCovMatrix}{\code{signature(object = "strumFittedModel")}: Accessor function, returns the fittedParametersCovMatrix of the object.}
    \item{chiTestOut}{\code{signature(object = "strumFittedModel")}: Accessor function, returns the chiTestOut of the object.}
    \item{show}{\code{signature(object = "strumFittedModel")}: Print a summary of the fitted model.}
  }
}
\seealso{
\code{\linkS4class{strumModel}}
}
\examples{
showClass("strumFittedModel")
}
\keyword{classes}
