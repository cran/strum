#==============================================================================
# File: computeDerivative.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver to calculate the derivatives of C and V at the optimum
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Estimate the 1st and 2nd derivatives of C and V
#------------------------------------------------------------------------------
.estimateDeltaDerivative = function(y, x, vc, betaMat, sigmaMat)
{
  numTrt = ncol(y$yAll[[1]])
  numCov = ncol(x$xAll[[1]])
  numVC  = length(vc$vcAll[[1]])
  numPed = length(y$yAll)
  
  fDeriv = matrix(ncol = numPed, nrow = 0)
  sDeriv = list()
  s1Deriv = list()
  s2Deriv = list()

  for( t1 in 1:numTrt )
  {	
    # 1. univariate derivatives
    #---------------------------

    bet1 = matrix(betaMat[,t1], ncol = 1)
    sig1 = lapply(sigmaMat, function(s) return(matrix(s[t1,t1])))

    # all
    #-----
    tmp_y1    = lapply(y$yAll, function(yk) return(yk[,t1]))
    missingY1 = lapply(tmp_y1, .findMissing)
    
    yt1  = mapply(.filterMissing,   tmp_y1,   missingY1, FALSE,     SIMPLIFY=FALSE) 
    xt1  = mapply(.filterMissing,   x$xAll,   missingY1, TRUE,      SIMPLIFY=FALSE) 
    vct1 = mapply(.filterMissingVC, vc$vcAll, missingY1, missingY1, SIMPLIFY=FALSE) 

    dUni = .Call("computeDeriv",
                 yt1, xt1, vct1, list(), list(), list(), list(), bet1, sig1)

    if( length(y$yPro) )
    {
      # probands
      #----------
      tmp_yp1    = lapply(y$yPro, function(yk) return(yk[,t1]))
      missingYp1 = lapply(tmp_yp1, .findMissing)
    
      ypt1  = mapply(.filterMissing,   tmp_yp1,  missingYp1, FALSE,      SIMPLIFY=FALSE) 
      xpt1  = mapply(.filterMissing,   x$xPro,   missingYp1, TRUE,       SIMPLIFY=FALSE) 
      vcpt1 = mapply(.filterMissingVC, vc$vcPro, missingYp1, missingYp1, SIMPLIFY=FALSE) 

      dUniP = .Call("computeDeriv",
                    ypt1, xpt1, vcpt1, list(), list(), list(), list(), bet1, sig1)

      fDeriv = rbind(fDeriv, dUni$dPar - dUniP$dPar)
      sDeriv[[length(sDeriv)+1]] = dUni$d2Par - dUniP$d2Par
    } else
    {
      fDeriv = rbind(fDeriv, dUni$dPar)
      sDeriv[[length(sDeriv)+1]] = dUni$d2Par
    }

    if( t1 == 1 )
      next

    # 2. bivariate derivatives
    #--------------------------
    for( t2 in (t1-1):1 )
    {
      bet12 = matrix(betaMat[,c(t1,t2)], ncol=2)
      sig12 = lapply(sigmaMat, function(s) return(s[c(t1,t2),c(t1,t2)]))

      # all
      #-----
      tmp_y2    = lapply(y$yAll, function(yk) return(yk[,t2]))
      missingY2 = lapply(tmp_y2, .findMissing)
    
      yt2  = mapply(.filterMissing,   tmp_y2,   missingY2, FALSE,     SIMPLIFY=FALSE)
      xt2  = mapply(.filterMissing,   x$xAll,   missingY2, TRUE,      SIMPLIFY=FALSE)
      vct2 = mapply(.filterMissingVC, vc$vcAll, missingY2, missingY2, SIMPLIFY=FALSE)
      vc12 = mapply(.filterMissingVC, vc$vcAll, missingY1, missingY2, SIMPLIFY=FALSE) 

      dBi = .Call("computeDeriv",
                  yt1, xt1, vct1, yt2, xt2, vct2, vc12, bet12, sig12)

      if( length(y$yPro) )
      {
        # probands
        #----------
        tmp_yp2    = lapply(y$yPro, function(yk) return(yk[,t2]))
        missingYp2 = lapply(tmp_yp2, .findMissing)
    
        ypt2  = mapply(.filterMissing,   tmp_yp2,   missingYp2, FALSE,      SIMPLIFY=FALSE)
        xpt2  = mapply(.filterMissing,   x$xPro,    missingYp2, TRUE,       SIMPLIFY=FALSE)
        vcpt2 = mapply(.filterMissingVC, vc$vcPro,  missingYp2, missingYp2, SIMPLIFY=FALSE)
        vcp12 = mapply(.filterMissingVC, vc$vcPro,  missingYp1, missingYp2, SIMPLIFY=FALSE) 

        dBiP = .Call("computeDeriv",
                     ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12, bet12, sig12)

        fDeriv = rbind(fDeriv, dBi$dPar - dBiP$dPar)
        sDeriv[[length(sDeriv)+1]]   = dBi$d2Par - dBiP$d2Par
        s1Deriv[[length(s1Deriv)+1]] = dBi$dP1dS - dBiP$dP1dS
        s2Deriv[[length(s2Deriv)+1]] = dBi$dP2dS - dBiP$dP2dS
      } else
      {
        fDeriv = rbind(fDeriv, dBi$dPar)
        sDeriv[[length(sDeriv)+1]]   = dBi$d2Par
        s1Deriv[[length(s1Deriv)+1]] = dBi$dP1dS
        s2Deriv[[length(s2Deriv)+1]] = dBi$dP2dS
      }      
    }
  }

  return(list(numTrt=numTrt, numCov=numCov, numVC=numVC,
              fDeriv=fDeriv, sDeriv=sDeriv, s1Deriv=s1Deriv, s2Deriv=s2Deriv))
}
