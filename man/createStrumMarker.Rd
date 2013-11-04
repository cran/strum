\name{createStrumMarker}
\alias{createStrumMarker}
\title{
 Create strumMarker Object
}
\description{
A function to create strumMarker object.
}
\usage{
createStrumMarker(hapMapData,
                  populationRecombRate = 50,
                  errorRate = 0,
                  mutationRate = 0,
                  missingRate = 0,
                  coding = c(0,1,2),
                  returnIBD = FALSE,
                  intervalIBD = 10)
}
\arguments{
  \item{hapMapData}{Object of class \code{data.frame} containing input data (see Details).}
  \item{populationRecombRate}{Numeric value stating the average number of 
  crossovers per cM when simulating founders.}
  \item{errorRate}{Numeric value stating the probability of a genotype error.  If 
  a genotype error occurs, a genotype (other than the true genotype) is selected 
  with uniform probability.}
  \item{mutationRate}{Numeric value stating the probability that an allele will 
  be mutated before it is passed on.  If a mutation occurs, a new allele (not 
  the original allele) is selected with uniform probability.}
  \item{missingRate}{Numeric value stating the probability of a missing genotype.}
  \item{coding}{Numeric vector of length three stating the coding of inheritance 
  model.  The first element corresponds to the code assigned to the homozygous 
  minor allele.  The second element corresponds to the heterozygous genotype.  
  The third element corresponds to the homozygous major allele.  The default 
  coding is additive (i.e. coding = c(0, 1, 2)).}
  \item{returnIBD}{Logical. If TRUE, then the ibd slot in the strumData object, 
  simulated using the returned strumMarker object from this function, will be 
  populated. If FALSE, the ibd slot will be left empty.}
  \item{intervalIBD}{Numeric value stating the interval that ibd values should be 
  calculated.}
}
\details{
The value of \code{hapMapData} must be a data.frame containing hapmap data.  The 
three columns that state the SNP information must present:
\enumerate{
  \item \code{rsID} contains the names of the marker.
  \item \code{chr} is the chromosome or linkage group.
  \item \code{phys_position} is the physical base pair position.  It is converted 
  to Haldane (genetic distance) map position in centiMorgans as phys_position/1000000.
}
}
\value{
Returns an object of class \code{strumMarker}.
}
\seealso{
\code{\linkS4class{strumMarker}}, \code{\link{importHapmapData}}
}
\examples{
\dontrun{
# Get some hapmap data & selct 10 snps.
#--------------------------------------
hap20 = importHapmapData(20)
hap20 = hap20[(1:10)*10,]

# Create strumMarker object with hapmap data.
#--------------------------------------------
snpStrumMarker = createStrumMarker(hap20)
}
}

