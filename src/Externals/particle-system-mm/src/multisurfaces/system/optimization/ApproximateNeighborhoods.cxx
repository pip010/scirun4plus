#include <cstdlib>
#include <iostream>
#include <list>
#include <math.h>
#include <time.h>
#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/domain/Domain.h>
#include <system/ParticleSystem.h>
#include <system/optimization/ApproximateNeighborhoods.h>

#ifdef _WIN32
#pragma warning( disable : 4244 4018 )
#endif

using namespace particle_sys;
using namespace custom_class;
using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
ApproximateNeighborhoods::ApproximateNeighborhoods( int intervals )
{
  _countdown = 0; // so that the first time around we build the lists
  _num_iterations = intervals;
}

//------------------------------------------------------------------------
// Function    : optimize
// Description : determine if the neighborhoods need to be built --
//               either because enough iterations have gone by, or
//               because some particles split or died -- and if so,
//               build the neighborhoods
//------------------------------------------------------------------------
void ApproximateNeighborhoods::optimize( svector<DynamicSurfacePoint*> &points )
{
  // check if we need to rebuild the neighborhoods
  if (1)// (_countdown == 0) ||
    //       (points[0]->domain()->neighborhood()->storageInt() == -1) )
  {
    //clock_t st;
    //st = clock();
    
    
    // set the flag in the neighborhood structure to indicate the 
    //   lists should be updated
    points[0]->domain()->neighborhood()->storageInt( 0 );

    //cout << "Rebuilding Neighborhoods [start]" << endl;

    // now, go through each point and build the neighborhood!
    for ( int i = 0; i < points.size(); i++ )
      points[i]->domain()->neighborhood()->determineNeighborhood(points[i]);

	//cout << "Rebuilding Neighborhoods [end]" << endl;

    // reset the countdown
    _countdown = _num_iterations-1;

    // and set the neighborhood to indicate the lists have been built
    points[0]->domain()->neighborhood()->storageInt( 1 );  
    //double seconds =  (double)(clock() - st)/(double)CLOCKS_PER_SEC;
    //cout << "ApproximateNeighborhoods::optimize iteration: "
    //     << seconds << " seconds." << endl;
  }

  // otherwise, don't update neighborhoods
  else
    --_countdown;

  _optimized = true;
}

