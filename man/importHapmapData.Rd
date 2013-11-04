\name{importHapmapData}
\alias{importHapmapData}
\title{
Import Hapmap Data
}
\description{
A function to import hapmap data from the website stated.
}
\usage{importHapmapData(chr, pop="CEU", ...)
}
\arguments{
  \item{chr}{Numeric value stating a chromosome number.}
  \item{pop}{Character value stating a population}
  \item{...}{Other arguments}
}
\source{
http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2009-02_phaseIII/HapMap3_r2/

CEU/TRIOS/hapmap3_r2_b36_fwd.consensus.qc.poly.\code{chr}_ceu.phased.gz
}
\examples{
\dontrun{
hap20 = importHapmapData(20)
}
}

