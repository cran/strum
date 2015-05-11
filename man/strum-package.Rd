\name{strum-package}
\docType{package}
\alias{strum-package}

\title{
A Package for STRUctural Modeling of Latent Variables for General Pedigree Data
}
\description{
This package implements the simulation and fitting of a broad range of latent 
variable and structural equation models.  It can handle multilevel models, polygenic 
random effects and linkage random effects.  Traditional structural equation models and 
confirmatory factor analysis may also be performed.  The framework implemented now can 
only handle quantitative variables.  Ordinal and nominal variables will be included 
in a future release.
}
\details{
\tabular{ll}{
Package: \tab strum\cr
Type: \tab Package\cr
Version: \tab 0.6.1\cr
Date: \tab 2015-05-11\cr
License: \tab GPL (>= 2)\cr
}
The following illustrates the typical steps for a strum analysis.

\enumerate{

\item Construct a strumModel by \code{\link{createStrumModel}} function.

\item Prepare data using \code{\link{createStrumData}} function.

\item Run analysis by the function call \code{\link{strum}}.
}

This package can also be used for simulation studies as following.
\enumerate{

\item Construct a strumSimModel by \code{\link{createSimModel}} function.

\item Simulate data using \code{\link{simulateStrumData}} function.

\item Run analysis on the simulated data by the function call \code{\link{strum}}.
}

Please check "Examples" section for the full illustration of strum use.
}
\examples{
\dontrun{
#================
# strum analysis
#================

# 1. Construct strumModel
#-------------------------

## 1.1 Model formulas
#---------------------
testForm1 = 'bp     =~ SBP + DBP
             anger  =~ A1 + A2
             stress =~ S1 + S2
             bp ~ anger + stress
             stress ~ anger + rs6040343   
             var(stress)=.1  
            '

testForm2 = 'L1 =~ SBP + DBP
             L1 ~ sex + <a,p,e>
            '

## 1.2 Create a strumModel
#--------------------------
myStrumModel1 = createStrumModel(formulas = testForm1)
myStrumModel2 = createStrumModel(formulas = testForm2)

# 2. Prepare data
#-----------------

## 2.1 Read a data file
#-----------------------
dF = read.table("simped.dat", header=T)

## 2.2 Create a strumData object
#--------------------------------

### 2.2.1 No IBD file
#---------------------
myStrumData1 = createStrumData(dF, "Pedigree")

### 2.2.2 With IBD file
#-----------------------
myStrumData2 = createStrumData(dF, "Pedigree", ibdFileName="GENIBD.chr1.ibd")

# 3. Run strum analysis
#-----------------------

## 3.1 Model with no ibd markers
#--------------------------------
myResult1 = strum(myStrumModel1, myStrumData1)

## 3.2 When an ibd marker is specified
#--------------------------------------
myResult2 = strum(myStrumModel2, myStrumData2, iMarkers=c("chr1marker1"))

#==================
# simulation study
#==================

# 1. Construct simModel
#-----------------------

## 1.1 Get some hapmap data & selct 10 snps
#-------------------------------------------
hap20 = importHapmapData(20) # 'load(file="hap20.rdata")' with saved hapmap data
hap20 = hap20[(1:10)*10,]

## 1.2 Create strumMarker object with hapmap data
#-------------------------------------------------
snpStrumMarker = createStrumMarker(hap20)

## 1.3 Ascertainment function
#-----------------------------
aFunction = function(data) return(any(data$disease == 1))

## 1.4 Model formula
#--------------------
simForm = 'bp     =~ SBP + DBP
           anger  =~ A1 + 0.5*A2
           stress =~ S1 + 0.9*S2
           bp ~ anger + stress + <p,e>
           stress ~ anger + rs6040343   
           var(stress)=.1  
          '

## 1.5 Create a strumModel
#--------------------------
mySimModel = createSimModel(
                 formulas = simForm,
                 markerInfo = snpStrumMarker,
                 ascertainment = aFunction
               )

# 2. Simulate data based on the given data structure 
#----------------------------------------------------

## 2.1 Read a data file
#-----------------------
dF = read.csv("chr1.csv")[,1:8]
names(dF) = c("family","id", "father","mother",names(dF)[5:8])
aStrumData = createStrumData(dF, "Pedigree")

## 2.2 Simulate data
#--------------------
mySimulatedStrumData = simulateStrumData(mySimModel, aStrumData)

# 3. Run strum analysis using simulated data
#--------------------------------------------
mySimResult = strum(myStrumModel1, mySimulatedStrumData)
}
}

\author{
Nathan Morris, Yeunjoo Song, Stephen Cahn

Maintainer: Nathan Morris <nathan.morris@cwru.edu>, Yeunjoo Song<yeunjoo.song@cwru.edu>
}
\references{
Morris, N.J., Elston, R.C., & Stein, C.M. (2010). A framework for structural equation 
models in general pedigrees. Human Heredity 70:278-286.

Song, Y.E., Stein, C.M., Morris, N.J. (2015). strum: an R package for \emph{stru}ctural \emph{m}odeling of latent variables for general pedigrees. BMC Genetics 16:35.
}
\keyword{ package }
