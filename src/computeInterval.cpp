#include "utils.h"

RcppExport SEXP computeInterval(SEXP numMarkers,
                                SEXP recomMapPos,
                                SEXP markerMapPos)
{
  colvec recomVec;

  try
  {
    int numM = as<int>(numMarkers);

    NumericVector rPos(recomMapPos);
    NumericVector mPos(markerMapPos);

    recomVec.zeros(rPos.size()*2);

    int currentStartPos = 0;
    int currentRecomPos = 0;

    for( int m = 0; m < numM; ++m )
    {
      if( mPos[m] >= rPos[currentRecomPos] )
      {
        recomVec[(currentRecomPos * 2)]     = currentStartPos+1;
        recomVec[(currentRecomPos * 2) + 1] = m;

        currentStartPos = m;
        currentRecomPos = currentRecomPos + 1;

        while( mPos[m] >= rPos[currentRecomPos] )
        {
          recomVec[(currentRecomPos * 2)]     = 0;
          recomVec[(currentRecomPos * 2) + 1] = 0;

          currentRecomPos = currentRecomPos + 1;
        }
      }
    }

    recomVec[(currentRecomPos * 2)]   = currentStartPos+1;
    recomVec[(currentRecomPos * 2)+1] = numM;
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return wrap(recomVec);
}
