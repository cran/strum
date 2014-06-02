#==============================================================================
# File: estimateParameter.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver to estimate the parameters.
#        - Calculate log-likelihood
#        - Calculate Delta parameters
#        - Calculate Theta parameters
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Calculate the model parameters in stage 1.
# - Estimate delta parameters, C and V in Equation 5 and 6.
#------------------------------------------------------------------------------
.estimateDeltaParameter = function(y, x, vc, step1OptimControl)
{
  numTrait = ncol(y$yAll[[1]])
  numCov   = ncol(x$xAll[[1]])
  numVC    = length(vc$vcAll[[1]])
  numPed   = length(y$yAll)

  betaMat  = matrix(0, nrow = numCov, ncol = numTrait, byrow = FALSE)
  sigmaMat = list()		
  for( j in 1:numVC )	
    sigmaMat[[j]] = matrix(0, numTrait, numTrait)
  
  ypt1 = list()
  ypt2 = list()
  xpt1 = list()
  xpt2 = list()
  vcpt1 = list()
  vcpt2 = list()
  vcp12 = list()

  for( t1 in 1:numTrait )
  {	
    #
    # 1. Compute univariate log-likelihood
    #--------------------------------------
    # all
    #-----
    tmp_y1    = lapply(y$yAll, function(yk) return(yk[,t1]))
    missingY1 = lapply(tmp_y1, .findMissing)
    
    yt1  = mapply(.filterMissing,   tmp_y1,   missingY1, FALSE,     SIMPLIFY=FALSE) 
    xt1  = mapply(.filterMissing,   x$xAll,   missingY1, TRUE,      SIMPLIFY=FALSE) 
    vct1 = mapply(.filterMissingVC, vc$vcAll, missingY1, missingY1, SIMPLIFY=FALSE) 

    if( length(y$yPro) )
    {
      # probands
      #----------
      tmp_yp1    = lapply(y$yPro, function(yk) return(yk[,t1]))
      missingYp1 = lapply(tmp_yp1, .findMissing)
    
      ypt1  = mapply(.filterMissing,   tmp_yp1,  missingYp1, FALSE,      SIMPLIFY=FALSE) 
      xpt1  = mapply(.filterMissing,   x$xPro,   missingYp1, TRUE,       SIMPLIFY=FALSE) 
      vcpt1 = mapply(.filterMissingVC, vc$vcPro, missingYp1, missingYp1, SIMPLIFY=FALSE) 
    }

    temp1 = .computeUnivariateLL(yt1, xt1, vct1, ypt1, xpt1, vcpt1, step1OptimControl)

    betaMat[,t1] = temp1$beta
    for( v in 1:numVC )
      sigmaMat[[v]][t1,t1] = temp1$sigma[v]

    if( t1 == 1 )
      next

    #
    # 2. Compute bivariate log-likelihood
    #-------------------------------------
    for( t2 in (t1-1):1 )
    {
      # all
      #-----
      tmp_y2    = lapply(y$yAll, function(yk) return(yk[,t2]))
      missingY2 = lapply(tmp_y2, .findMissing)
    
      yt2  = mapply(.filterMissing,   tmp_y2,   missingY2, FALSE,     SIMPLIFY=FALSE)
      xt2  = mapply(.filterMissing,   x$xAll,   missingY2, TRUE,      SIMPLIFY=FALSE)
      vct2 = mapply(.filterMissingVC, vc$vcAll, missingY2, missingY2, SIMPLIFY=FALSE)
      vc12 = mapply(.filterMissingVC, vc$vcAll, missingY1, missingY2, SIMPLIFY=FALSE) 

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
      }

      bet12 = matrix(betaMat[,c(t1,t2)], ncol=2)
      sig12 = lapply(sigmaMat, function(s) return(diag(s)[c(t1,t2)]))

      temp2 = .computeBivariateLL(yt1,  xt1,  vct1,  yt2,  xt2,  vct2,  vc12,
                                 ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12,
                                 bet12, sig12, step1OptimControl)

      for( v in 1:numVC )	
      {	
        sigmaMat[[v]][t1,t2] = temp2$sigma[v]
        sigmaMat[[v]][t2,t1] = temp2$sigma[v]
      }
    }
  }

  return(list(C=betaMat, V=sigmaMat))
}

#------------------------------------------------------------------------------
# Calculate the model parameters in stage 2.
# - Estimate reduced model parameters (theta), Equation 10.
#------------------------------------------------------------------------------
.estimateThetaParameter = function(model, delta, w, startTheta, step2OptimControl)
{
  positive = grep("<", paramNames(model))

  lengthNV = length(paramNames(model)) - length(positive)

  Qstar = function(thetaNV)
          {
            mMat = .getModelMats(model, thetaNV)

            B  = mMat$B
            L  = mMat$L
            Gs = mMat$Gs
            Gm = mMat$Gm

            Zs = mMat$Zs
            Es = mMat$Es

            q = .Call("computeQstar", B, L, Gs, Gm, Zs, Es, as.matrix(delta), as.matrix(w))

            return(q)
          }

  optimOut = list()
  mymin = 10^100

  for( i in 1:length(startTheta) )
  {
    optimOut[[i]] = tryCatch(
                      optim(startTheta[[i]]$startVal[1:lengthNV],
                            Qstar,
                            control = step2OptimControl,
                            method = "BFGS"),
                      error = function(d)
                              {
                                sval = startTheta[[i]]$startVal[1:lengthNV] + rnorm(lengthNV, sd=.1)
                                return(tryCatch(
                                         optim(sval,
                                               Qstar,
                                               control = step2OptimControl,
                                               method = "BFGS"),
                                         error = function(d) list(value=10^100,par=NA)))
                              })

    if( abs((optimOut[[i]]$value-mymin)/mymin) < 10^-4 )
      break

    if( optimOut[[i]]$value < mymin )
      mymin = optimOut[[i]]$value
  }
  
  tmp = sapply(optimOut, function(ooi) ooi$value)
  minOut = c(which(tmp==min(tmp)))
  optimOut = optimOut[[minOut[1]]]

  if( optimOut$convergence != 0 )
    warning("optim did not converge in stage 2.")

  mMat = .getModelMats(model, optimOut$par)

  B  = mMat$B
  L  = mMat$L
  Gs = mMat$Gs
  Gm = mMat$Gm

  Zs = mMat$Zs
  Es = mMat$Es

  thetaV = .Call("computeThetaV", B, L, Gs, Gm, Zs, Es, as.matrix(delta), as.matrix(w))

  optimPar = c(optimOut$par, thetaV)
  names(optimPar) = paramNames(model)

  return(optimPar)
}

#------------------------------------------------------------------------------
# A function which utilizes the C likelihood function to estimate
#  the log-likelihood of "saturated model" - univariate .
#------------------------------------------------------------------------------
.computeUnivariateLL = function(y, x, vc, yp, xp, vcp, step1OptimControl)
{
  numCov = ncol(x[[1]])
  numVC  = length(vc[[1]])         
  
  yb = do.call("c",y)
  xb = do.call(rbind,x)

  fit = lm(yb ~ xb-1)
  betaStart = coef(fit)
  sigmaTot = summary(fit)$sigma

  startVal = c(rep(sigmaTot,numVC)/numVC, betaStart)

  logLikelihoodP = 0.0;

  Q = function(theta)
      {
        sigma = .makeSigmaList(theta[1:numVC], list())

        pd = sapply(sigma, .testPosDefMatrix)
        if( any(pd == FALSE) )
          return(-1*10^6)

        beta = matrix(theta[(numVC+1):(numVC + numCov)], ncol = 1)

        logLikelihood = .Call("computeLL", 
                              y, x, vc, list(), list(), list(), list(), beta, sigma)

        if( length(yp) )
          logLikelihoodP = .Call("computeLL", 
                                 yp, xp, vcp, list(), list(), list(), list(), beta, sigma)

        return(logLikelihood - logLikelihoodP)
      }
  
  optimOut = optim(startVal, Q, control=step1OptimControl)

  if( optimOut$convergence != 0 )
    warning("optim did not converge in stage 1.")
  
  ret = list(sigma = optimOut$par[1:numVC],
             beta  = optimOut$par[(numVC+1):(numVC + numCov)])

  return(ret)	
}

#------------------------------------------------------------------------------
# A function which utilizes the C likelihood function to estimate
#  the log-likelihood of "saturated model" - bivariate .
#------------------------------------------------------------------------------
.computeBivariateLL = function(yt1,  xt1,  vct1,  yt2,  xt2,  vct2,  vc12,
                               ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12,
                               beta, sig, step1OptimControl)
{
  numVC = length(vct1[[1]])         

  startVal1 = rep(0, numVC)  
  startVal2 = .rdirichlet(1, rep(1, numVC))

  logLikelihoodP = 0.0;

  Q = function(theta)
      {
        sigma = .makeSigmaList(sig, theta)

        pd = sapply(sigma, .testPosDefMatrix)
        if( any(pd == FALSE) )
          return(-1*10^6)

        logLikelihood = .Call("computeLL", 
                              yt1, xt1, vct1, yt2, xt2, vct2, vc12, beta, sigma)

        if( length(ypt1) )
          logLikelihoodP = .Call("computeLL", 
                                 ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12, beta, sigma)

        return(logLikelihood - logLikelihoodP)
      }

  optimOut1 = optim(startVal1, Q, control=step1OptimControl)

  if( optimOut1$convergence != 0 )
    warning("optim did not converge in stage 1.")

  optimOut2 = optim(startVal2, Q, control=step1OptimControl)

  if( optimOut2$convergence != 0 )
    warning("optim did not converge in stage 1.")

  logLLs   = c(optimOut1$value, optimOut2$value)
  maxLLPos = which.max(logLLs)

  optimPars = rbind(optimOut1$par, optimOut2$par)

  sigVal0 = sqrt(sig[[1]][1] * sig[[1]][2]) * optimPars[maxLLPos,][1]
  if( numVC > 1 )
    for( v in 2:numVC )
      sigVal0 = c(sigVal0, sqrt(sig[[v]][1] * sig[[v]][2]) * optimPars[maxLLPos,][v])

  ret = list(sigma = sigVal0)

  return(ret)	
}

#------------------------------------------------------------------------------
# Check whether a matrix is positive definite or not.
#------------------------------------------------------------------------------
.testPosDefMatrix = function(mat)
{
  if( nrow(mat) == 1 )
  {
    if( mat[1,1] > 0 )
      return(TRUE)
    else
      return(FALSE)
  }

  if( (mat[1,1]>0) & (mat[2,2]>0) & ((mat[1,1]*mat[2,2]) > mat[2,1]*mat[1,2]) )
    return(TRUE)
  else  
    return(FALSE)
}

#------------------------------------------------------------------------------
# Form a list of variacne-covariance matrices.
#------------------------------------------------------------------------------
.makeSigmaList = function(vari, cor)
{
  numVC = length(vari)
  sigmaL = list()

  if( length(cor) == 0 )
    for( v in 1:numVC )   
      sigmaL[[v]] = matrix(c(vari[v]), nrow = 1, ncol = 1, byrow = FALSE)
  else
    for( v in 1:numVC )
    {
        covari = sqrt(vari[[v]][1]*vari[[v]][2]) * cor[v]

        sigmaL[[v]] = matrix(c(vari[[v]][1],covari,covari,vari[[v]][2]),
                             nrow = 2, ncol = 2, byrow = FALSE)
    }

  return(sigmaL)
}

#------------------------------------------------------------------------------
# Get the model matrices
#------------------------------------------------------------------------------
.getModelMats = function(model, thetaNV)
{
  lengthNV = length(thetaNV)
  lengthV  = length(paramNames(model)) - lengthNV
  ei= diag(1, lengthV)

  thBase = c(thetaNV, rep(0,lengthV))

  B  = model@B(thBase)
  L  = model@L(thBase)
  Gs = model@Gs(thBase)
  Gm = model@Gm(thBase)

  Zs = list()
  Es = list()

  Z = lapply(model@Z, function(v) v(thBase))
  E = lapply(model@E, function(v) v(thBase))

  Zs[[1]] = Z
  Es[[1]] = E

  for( i in 1:lengthV )
  {
    thAll = c(thetaNV, ei[i,])

    Z = lapply(model@Z, function(v) v(thAll))
    E = lapply(model@E, function(v) v(thAll))

    Zs[[i+1]] = Z
    Es[[i+1]] = E
  }

  return(list(B=B, L=L, Gs=Gs, Gm=Gm, Zs=Zs, Es=Es))
}
