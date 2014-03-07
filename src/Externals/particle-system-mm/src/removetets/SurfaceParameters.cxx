#include <iostream>
#include <cstdlib>
#include "SurfaceParameters.h"

using namespace std;

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the SurfaceParams
//------------------------------------------------------------------------
ostream& operator << ( ostream& outs, const SurfacePointParams& source )
{
  outs << "              " << 
    "<F,Fx,Fxx> : <" << source._F <<
    "," << source._Fx << "\n" << source._Fxx;
  
  return outs; 
};

