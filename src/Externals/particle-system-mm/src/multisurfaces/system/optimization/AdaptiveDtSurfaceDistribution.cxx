#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/domain/Domain.h>
#include <system/ParticleSystem.h>
#include <system/optimization/SurfacePopulationController.h>
#include <system/optimization/AdaptiveDtSurfaceDistribution.h>

//C++11
#include <thread>
#include <functional>
#include <mutex>

#ifdef _WIN32
#pragma warning (disable : 4244)
#endif

using namespace particle_sys;
using namespace custom_class;
using namespace std;

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : 
//------------------------------------------------------------------------
AdaptiveDtSurfaceDistribution::AdaptiveDtSurfaceDistribution( float 
                                                              F_threshold )
{
  _constraint = new SurfaceConstraint;

  _prev_global_energy = 0.0;
  _num_iterations = 0;
  _F_threshold = F_threshold;
}

AdaptiveDtSurfaceDistribution::~AdaptiveDtSurfaceDistribution()
{
  delete _constraint;
}

void parallelKernel(std::vector<DynamicSurfacePoint*> &points,SurfaceConstraint* constraint )
{
	for ( unsigned j = 0; j < points.size(); j++ )
	{
	  constraint->projectOntoSurface( points[j] );
	  
	  if ( points[j]->moved_outside() )
	  {
		 cout << "[[[SHABANG]]]" ;
		(points[j]->system())->removePoint( j );
		  j--;
	  }
	}
} 

//------------------------------------------------------------------------
// Function    : init
// Description : initialize the Points by projecting them onto the Surface
//------------------------------------------------------------------------
void AdaptiveDtSurfaceDistribution::init( 
  svector<DynamicSurfacePoint*> &points, int num_iterations )
{
  //P.Petrov 2014 attempt at parallelization
  //Thread t1;
  //Thread t2;
/*
  unsigned int  size = points.size();
  std::vector<DynamicSurfacePoint*> t1v(points._vector.begin(), points._vector.begin() + points._vector.size()/2);
  std::vector<DynamicSurfacePoint*> t2v(points._vector.begin() + points._vector.size()/2, points._vector.end());
  
  while ( num_iterations-- )
  {

  
	if(size > 1000)
    {
	  cout << "Recommend to go parallel for: " << t1v.size() <<"  "<< t2v.size() << "|||" << num_iterations << endl;
	
	
	  thread t1(parallelKernel,std::ref(t1v),_constraint);
	  thread t2(parallelKernel,std::ref(t2v),_constraint);
	  t1.join();
	  t2.join();
	  
	  cout << "\nNEXT ITER\n";
	}
	//else if( size > 100)
	//{
	//  thread t1(parallelKernel,std::ref(points),_constraint,0,size);
	//  t1.join();
	//}
	else
	{
		parallelKernel(points._vector,_constraint);
	}
	
	for ( unsigned j = 0; j < size; j++ )
	{
	  // check if the Point needs to be removed from the System
	  if ( points[j]->moved_outside() )
	  {
		(points[j]->system())->removePoint( j );
		  j--;
	  }
	}
	
 }
    */
    

   //OLD CODE
  while ( num_iterations-- )
  {
    for ( unsigned j = 0; j < points.size(); j++ )
    {
		//cout << " I"; 
      _constraint->projectOntoSurface( points[j] );
		
      // check if the Point needs to be removed from the System
      if ( points[j]->moved_outside() )
      {
        (points[j]->system())->removePoint( j );
	      j--;
	      //cout<< " R";
      }
    }
}
 
 /*
 cout << endl;
 
     for ( unsigned j = 0; j < points.size(); j++ )
    {
		cout << *(points[j]) << endl;
	}
  
 cout << endl;
 */
  // now, make sure the points are all within the threshold of the
  //   surface
  if ( !points.empty() )
  {
    //cout << points.size() << endl;
    points[0]->system()->cleanUpSystem( _F_threshold );
  }
}

//------------------------------------------------------------------------
// Function    : optimize
// Description : 
//------------------------------------------------------------------------
void 
AdaptiveDtSurfaceDistribution::optimize(svector<DynamicSurfacePoint*> &points)
{
  float energy, lambda, max_movement, energy_tmp;
  vector_type force, dx, old_position, old_normal;
  bool keep_iterating;
  clock_t st, elapsed, tmp;
  st = clock();
  elapsed = st;

  double s1, s2, s3; //3 sections to profile inside iteration loop
  s1 = s2 = s3 = 0.0l; //3 sections to profile inside iteration loop
  size_t size_pnts = points.size();

  vector<int> delP;
  for (unsigned i = 0; i < points.size(); ++i)
  {
    // get the original Point values
    lambda = points[i]->lambda();
    old_normal = points[i]->normal();
    old_position = points[i]->position();

    // get the energy and force at this Point
    points[i]->computeEnergyForce( energy, force );

    // compute Point movement in the local tangent plane
    dx = (-force) / lambda;
    //dx -= DotProduct(old_normal, dx) * old_normal;  
    points[i]->domain()->surface()->projectMotion( points[i]->position(),
                                                   old_normal, dx );

    // profiling the function...
    tmp = clock();
    s1 += (double)(tmp - elapsed)/(double)CLOCKS_PER_SEC;
    elapsed = tmp;

    // check if the movement is too big
    max_movement = 2.0*points[i]->radius();
    if ( dx.length() > max_movement) {
      energy_tmp = MAX_VALUE;
      //     cerr << "not moving: dx.len:" << dx.length() << ", max:" 
      //    << max_movement << endl;
      // otherwise, move the Point
    } else {
      points[i]->move( old_position+dx );

      if ( points[i]->moved_outside() )
        energy_tmp = MAX_VALUE;
      else
      {
        _constraint->projectOntoSurface(points[i], _F_threshold);

        // check if it moved outside the domain
        if ( points[i]->moved_outside() )
          energy_tmp = MAX_VALUE;

        // check if the particle has moved off of the surface
        else if ( fabs(points[i]->F()) > _F_threshold )
          energy_tmp = MAX_VALUE;

        // otherwise, lets get the new energy!
        else
          points[i]->computeEnergy( energy_tmp );
      }
    }

    //
    // TWO CASES
    //
    //   1 - energy not better with move, so start iterating and 
    //       increasing lambda until we get to a good value -- this
    //       means that we are taking too big of a step and missing
    //       the minimum, so make the step smaller and smaller until we
    //       no longer miss the min
    //
    //   2 - the energy is smaller with this move, so take the move and
    //       and decrease lambda -- this means that we are moving in
    //       the right direction, so lets try a little bigger step
    //       size next time
    //    

    // profiling the function...
    tmp = clock();
    s2 += (double)(tmp - elapsed)/(double)CLOCKS_PER_SEC;
    elapsed = tmp;

    // CASE 1
    if ( energy_tmp > energy ) 
    {    
      // first check to see if the lambda is already big, if so, just
      //   stay where it is and don't do anything!
      if ( lambda >= MAX_LAMBDA )
        // reset the position
        points[i]->move( old_position );

      // otherwise, lets keep making lambda bigger until we get to
      //   a good energy
      else
      {
        keep_iterating = true;
        while ( keep_iterating )
        {
          energy_tmp = 0.0;

          // increase lambda
          lambda *= 10.0;

          // get the new (smaller) dx in the tangent plane
          dx = (-force) / lambda;
          //dx -= DotProduct(old_normal, dx) * old_normal;
          points[i]->domain()->surface()->projectMotion(points[i]->position(),
                                                        old_normal, dx);

          // check if the movement is too big
          if ( dx.length() > max_movement )
            energy_tmp = MAX_VALUE;

          // otherwise, move the Point
          else
          {
            points[i]->move( old_position+dx );

            if ( points[i]->moved_outside() )
              energy_tmp = MAX_VALUE;
            else
            {
              _constraint->projectOntoSurface( points[i] );

              // check if it moved outside the domain
              if ( points[i]->moved_outside() )
                energy_tmp = MAX_VALUE;

              //// check if the particle has moved off of the surface
              //else if ( fabs(points[i]->F()) > _F_threshold )
              //  energy_tmp = MAX_VALUE;

              // otherwise, lets get the new energy!
              else
                points[i]->computeEnergy( energy_tmp );
            }
          }
          
          // if the energy is finally lowered and the particle is inside 
          //   the domain, then lets take the move and quit
          if ( (energy_tmp < energy) && !points[i]->moved_outside() )
          {
            points[i]->lambda( lambda );
            keep_iterating = false;
          }

          // however, if the lambda is huge we should just keep the point
          //    at its original position (keeps particles on boundaries
          //    from constantly dying)
          else if ( lambda >= MAX_LAMBDA )
          {
            points[i]->move( old_position );
            points[i]->lambda( lambda );
            keep_iterating = false;

            // this is for some sort of weirdness with the fe stuff!
            if (points[i]->moved_outside()) {
              delP.push_back(i);
            }
          }
      
          // otherwise, let's decrease the step size and try again!

        } // while ( keep_iterating ) ...
      }
    } // if ( energy_tmp > energy ) ...

    // CASE 2
    else 
    {
      // do a check for lambda ~= 0 to keep things from infinitly
      //   equaling zero
      if ( lambda > MIN_LAMBDA )
        lambda *= 0.1;

      points[i]->lambda( lambda );
    }

    if (fabs(points[i]->F()) > _F_threshold) {
      delP.push_back(i);
    }
    // profiling the function...
    tmp = clock();
    s3 += (double)(tmp - elapsed)/(double)CLOCKS_PER_SEC;
    elapsed = tmp;

  } // for ( int i = 0; i < points.size(); i++ ) ...

  // delete the points flagged for deletion.
  sort(delP.begin(), delP.end());
  vector<int>::reverse_iterator iter = delP.rbegin();
  while (iter != delP.rend()) {
    int i = *iter++;
    (points[i]->system())->removePoint(i);
  }

  _optimized = doneMoving( points );
  ++_num_iterations;

  //double seconds =  (double)(clock() - st)/(double)CLOCKS_PER_SEC;

  //cout << "section 1: "<< s1 / seconds * 100. << "%" << endl;
  //cout << "section 2: "<< s2 / seconds * 100. << "%" << endl;
  //cout << "section 3: "<< s3 / seconds * 100. << "%" << endl;


  //cout << "AdaptiveDtSurfaceDistribution::optimize iteration: "
  //     << seconds << " seconds." << endl;
}

//------------------------------------------------------------------------
// Function    : doneMoving
// Description : check whether the system is at an steady state, return
//               true if so... --> this is based on the change of system
//               energy
//------------------------------------------------------------------------
bool AdaptiveDtSurfaceDistribution::doneMoving( svector<DynamicSurfacePoint*> 
                                                &points )
{
  if ( _num_iterations > 50 )
  {
    _num_iterations = 0;
    return true;
  }

  // get the new global energy of the Points
  float global_energy = 0.0;
  for ( unsigned i = 0; i < points.size(); i++ )
    global_energy += points[i]->energy();

  float change = fabs(global_energy - _prev_global_energy) / 
                 (_prev_global_energy + EPSILON);
  _prev_global_energy = global_energy;

  // now, check to see if the energy has changed much
  if ( change > 0.00015 ){
    return false;
  }
  // set the lambda values of the Points for the next time around
  for ( unsigned i = 0; i < points.size(); i++ )
    points[i]->lambda( INITIAL_LAMBDA );

  _num_iterations = 0;
  return true;
}




  

