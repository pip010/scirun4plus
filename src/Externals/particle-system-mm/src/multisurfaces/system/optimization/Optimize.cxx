#include <cstdlib>
#include <iostream>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/optimization/Optimize.h>
#include <ctime>

using namespace particle_sys;
using namespace std;

void tprint()
{
  time_t rawtime;
  //struct tm * ptm;

  time ( &rawtime );

  //ptm = gmtime ( &rawtime );

cout << "DEBUG TIMESTAMP " << rawtime << endl;

}
//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
Optimize::Optimize( Optimization **ops, int num_ops )
{
  // the number of optimizations the System has
  _num_ops = num_ops;

  // now, copy the location of the array of Optimizations
  _optimizations = ops;
}

Optimize::~Optimize()
{
  for ( int i = 0; i < _num_ops; i++ )
    delete _optimizations[i];
  delete [] _optimizations;
}

//------------------------------------------------------------------------
// Function    : init
// Description : 
//------------------------------------------------------------------------
void Optimize::init( custom_class::svector<DynamicSurfacePoint*> &points,
                    int num_iterations )
{
  for ( int i = 0; i < _num_ops; i++ )
    _optimizations[i]->init( points, num_iterations );
}

//------------------------------------------------------------------------
// Function    : optimize
// Description : this function does one step of optimization of the 
//               System --> the list of Optimizations is gone through
//               in order, with subsequent Optimizations only being
//               called if the previous one was optimized
//------------------------------------------------------------------------
void Optimize::optimize( svector<DynamicSurfacePoint*> &points )
{
  // first, check if we have no optimizations to do
  if ( !_num_ops )
    return;

  // otherwise, lets set the optimized flag to false, and see if gets 
  //   changed
  _optimized = false;

  for ( int i = 0; i < _num_ops; i++ )
  {
	  //DEBUG pip
	  //cout << "DEBUG START op " << i << endl;
	  //tprint();
	  
    _optimizations[i]->optimize( points );
    
	//DEBUG pip
	  //cout << "DEBUG END op " << i  << endl;
	  //tprint();
    
    if ( !_optimizations[i]->optimized()  )
      return;
  }

  // seems that all of the optimizations are finished, so set the 
  //   optimized flag back to true
  _optimized = true;
}




  

