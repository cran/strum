\name{strumModel-class}
\docType{class}

\alias{strumModel-class}
\alias{show,strumModel-method}

\title{Class \code{"strumModel"}}
\description{
\code{strumModel} is an S4 class that represents a STRUctural Modeling of latent 
Variables for family data.
}
\section{Objects from the Class}{
Objects should not be created by calls of the form \code{new("strumModel", ...)} 
but by the calls to the \code{\link{createStrumModel}} function. 
}
\section{Slots}{
  \describe{
    \item{\code{varList}:}{Object of class \code{"data.frame"} containing a list of variables and their properties.}
    \item{\code{formulas}:}{Object of class \code{"character"} containing the information about the formulas and dependencies.}
    \item{\code{allRandomEffects}:}{Object of class \code{"character"} listing all variance components in the model.}
    \item{\code{paramNames}:}{Object of class \code{"character"} describing the names of the model parameters.}
    \item{\code{ascertainment}:}{Object of class \code{"ANY"} stating the ascertainment criteria (function) of the data.}
    \item{\code{E}:}{Object of class \code{"list"} containing the variance matrices of measurement errors of the model.}
    \item{\code{Z}:}{Object of class \code{"list"} containing the variance matrices of random effects for the structural equation of the model.}
    \item{\code{L}:}{Object of class \code{"function"} to generate a matrix to relate the unobserved latent factors to the observed traits in the model.}
    \item{\code{B}:}{Object of class \code{"function"} to generate a matrix denoting the causal relationship among the latent factors in the model.}
    \item{\code{Gs}:}{Object of class \code{"function"} to generate a matrix to relate the observed covariates to the unobserved latent factors in the model.}
    \item{\code{Gm}:}{Object of class \code{"function"} to generate a matrix to relate the observed covariates to the observed traits in the model.}
    \item{\code{thToThB}:}{Object of class \code{"numeric"} for future use.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{strumVirtualModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "strumModel")}: Print a summary of the model.}
  }
}
\seealso{
\code{\linkS4class{strumSimModel}},  \code{\linkS4class{strumVirtualModel}},
\code{\link{createStrumModel}}
}
\examples{
showClass("strumModel")
}
\keyword{classes}
