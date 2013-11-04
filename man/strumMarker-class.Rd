\name{strumMarker-class}
\docType{class}

\alias{strumMarker-class}
\alias{coding<-}
\alias{coding}
\alias{errorRate<-}
\alias{errorRate}
\alias{haplotypes<-}
\alias{haplotypes}
\alias{intervalIBD<-}
\alias{intervalIBD}
\alias{markerFacts<-}
\alias{markerFacts}
\alias{missingRate<-}
\alias{missingRate}
\alias{mutationRate<-}
\alias{mutationRate}
\alias{populationRecombRate<-}
\alias{populationRecombRate}
\alias{returnIBD<-}
\alias{returnIBD}
\alias{coding<-,strumMarker-method}
\alias{coding,strumMarker-method}
\alias{errorRate<-,strumMarker-method}
\alias{errorRate,strumMarker-method}
\alias{haplotypes<-,strumMarker-method}
\alias{haplotypes,strumMarker-method}
\alias{intervalIBD<-,strumMarker-method}
\alias{intervalIBD,strumMarker-method}
\alias{markerFacts<-,strumMarker-method}
\alias{markerFacts,strumMarker-method}
\alias{missingRate<-,strumMarker-method}
\alias{missingRate,strumMarker-method}
\alias{mutationRate<-,strumMarker-method}
\alias{mutationRate,strumMarker-method}
\alias{populationRecombRate<-,strumMarker-method}
\alias{populationRecombRate,strumMarker-method}
\alias{returnIBD<-,strumMarker-method}
\alias{returnIBD,strumMarker-method}
\alias{show,strumMarker-method}

\title{Class \code{"strumMarker"}}
\description{
\code{strumMarker} is an S4 class that represents a marker data for a STRUctural Modeling of latent variables for family data.
}
\section{Objects from the Class}{
Objects should not be created by calls of the form \code{new("strumMarker", ...)} but by the calls to the \code{\link{createStrumMarker}} function. 
}

\section{Slots}{
  \describe{
    \item{\code{markerFacts}:}{Object of class \code{"data.frame"} containing the SNP information - markerName, chrom, mapPos, minorAllele, majorAllele.}
    \item{\code{haplotypes}:}{Object of class \code{"matrix"} containing hapmap data.}
    \item{\code{populationRecombRate}:}{Numeric value stating the average number of crossovers per cM when simulating founders.}
    \item{\code{errorRate}:}{Numeric value stating the probability of a genotype error.  If a genotype error occurs, a genotype (other than the true genotype) is selected with uniform probability.}
    \item{\code{mutationRate}:}{Numeric value stating the probability that an allele will be mutated before it is passed on.  If a mutation occurs, a new allele (not the original allele) is selected with uniform probability.}
    \item{\code{missingRate}:}{Numeric value stating the probability of a missing genotype.}
    \item{\code{coding}:}{Numeric vector of length three stating the coding of inheritance model.  The first element corresponds to the code assigned to the homozygous minor allele.  The second element corresponds to the heterozygous genotype.  The third element corresponds to the homozygous major allele.  The default coding is additive (i.e. coding = c(0, 1, 2)).}
    \item{\code{returnIBD}:}{Logical. If TRUE, then the ibd slot in the strumData object, simulated using the returned strumMarker object from this function, will be populated.  If FALSE, the ibd slot will be left empty.}
    \item{\code{intervalIBD}:}{Numeric value stating the interval that ibd values should be calculated.}
  }
}
\section{Methods}{
  \describe{
    \item{markerFacts}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{markerFacts}.}
    \item{haplotypes}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{haplotypes}.}
    \item{populationRecombRate}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{populationRecombRate}.}
    \item{errorRate}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{errorRate}.}
    \item{mutationRate}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{mutationRate}.}
    \item{missingRate}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{missingRate}.}
    \item{coding}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{coding}.}
    \item{returnIBD}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{returnIBD}.}
    \item{intervalIBD}{\code{signature(object = "strumMarker")}: Accessor function, returns the value of the slot \code{intervalIBD}.}
    \item{show}{\code{signature(object = "strumMarker")}: Print a summary of the class.}
  }
}
\seealso{
\code{\link{createStrumMarker}},  \code{\link{simulateStrumData}}
}
\examples{
showClass("strumMarker")
}
\keyword{classes}
