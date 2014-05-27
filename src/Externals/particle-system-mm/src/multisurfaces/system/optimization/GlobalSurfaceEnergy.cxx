#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <list>
#include <math.h>
#include <time.h>
#include <features/mtxlib.h>
#include <system/defines.h>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/domain/Domain.h>
#include <system/ParticleSystem.h>
#include <system/optimization/GlobalSurfaceEnergy.h>

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
GlobalSurfaceEnergy::GlobalSurfaceEnergy(bool biased_splitting_dying,
                                         float threshold) :
  Optimization(), SurfacePopulationController(threshold)
{
  _constraint = new SurfaceConstraint;
  _biased_splitting_dying = biased_splitting_dying;
  _first_time_through = true;
}

GlobalSurfaceEnergy::~GlobalSurfaceEnergy()
{
  delete _constraint;
}

//------------------------------------------------------------------------
// Function    : optimize
// Description : 
//------------------------------------------------------------------------
void GlobalSurfaceEnergy::optimize( svector<DynamicSurfacePoint*> &points )
{
  cout << "GlobalSurfaceEnergyNA optmize" << endl;
	
  if ( _biased_splitting_dying )
    biasedOptimize( points );
  else
    randomOptimize( points );

  //_optimized = true;
}

/************************************************************************/

//------------------------------------------------------------------------
// Function    : biasedOptimize
// Description : optimize the system using a biased splitting and dying
//               algorithm
//------------------------------------------------------------------------
void GlobalSurfaceEnergyNA::biasedOptimize( svector<DynamicSurfacePoint*> 
                                            &points )
{
  clock_t st;
  st = clock();
  _optimized = true;

  bool force_op;
  if ( _first_time_through ) force_op = true;
  else force_op = false;
  _first_time_through = false;

  if ( points.empty() )
    return;

  

  if ( points[0]->domain()->numNeighbors() == 0 )
  {
    int num_pts = points.size();
    float local_energy;
    for ( int i = 0; i < num_pts; i++ )
    {
      points[i]->computeEnergy( local_energy );

      if ( local_energy != 0.0 )
      {
        deleteAPoint( i, points[i] );
        --i;
        --num_pts;
      }
    }

    return;
  }

  // do a first pass to compute all the local energies -- seems to work
  //   better if we compute them all initially then inside of the
  //   splitting and deleting loop!
  float local_energy, bias;
  int num_split = 0, num_deleted = 0;

  // check for quad points which should be singular!!!
  float ideal_local_energy = points[0]->domain()->idealEnergy();

  float delete_val, split_val;
  if ( points[0]->domain()->numNeighbors() == 2 )
  {
    delete_val = 2.5;
    split_val = -0.25;
  }
  else
  {
    delete_val = 0.75;
    split_val = -0.35;
  }

  for ( int i = 0; i < points.size(); i++ )
  {
    points[i]->computeEnergy( local_energy );

    // bias is the percentage of the ideal energy that the Particle
    //   is missing or has extra
    if ( points[i]->energy() != 0.0 )
      bias = (rand()/(float)RAND_MAX) * 
        (local_energy - ideal_local_energy) / ideal_local_energy; 
    else
      bias = -1.0;

    // store it in a temporary variable
    points[i]->storageFloat( bias );

    // count the number of particles that will be split or deleted
    if ( (bias < split_val) || (bias > delete_val) )
    {

      // let's see how much we change the population
      if ( bias > 0.0 )
        ++num_deleted;
      else
        ++num_split;
    }
  }

  // if less than 1% of the particles will be split or deleted,
  //   then don't do it!
  if ( points[0]->domain()->numNeighbors() == 6 )
  {
    if ( (((float)(num_split+num_deleted)/(float)points.size()) < 0.01) &&
         !force_op )
    //if ( (((float)num_split/(float)points.size()) < 0.01) && 
    //     (((float)num_deleted/(float)points.size()) < 0.01) )
    {
      cout << num_split << "  " << num_deleted << endl;
      return;
    }

    //if ( (fabs((double)(num_split-num_deleted))/
    //    (float)(num_split+num_deleted)) < 0.3 )
    //  return;
  }
  else if ( points[0]->domain()->numNeighbors() == 2 )
  {
    if ( (((float)(num_split+num_deleted)/(float)points.size()) < 0.005) &&
         !force_op )
      return;

     
    if ( ((fabs((double)(num_split-num_deleted))/
        (float)(num_split+num_deleted)) < 0.1) && !force_op )
      return;
  }

  // go through the Points, and do a biased splitting and deleting
  //   based on the local energy value at the Points

  cout << endl;
  cout << "Num split = " << num_split << "    Num delete = " << 
    num_deleted << "  " << points.size() << endl << endl;

  vector<int> delP;
  vector<int> s4WSF;
  vector<int> s4;
  vector<int> split;

  int num_pts = points.size();
  num_split=0; num_deleted=0;
  for ( int i = 0; i < num_pts; i++ )
  {  
    //
    // REMEMBER --> bias is based on the difference between the
    //              local energy and the ideal energy; a small
    //              bias means that the density of the neighborhood
    //              is close to ideal!!!
    //
    bias = points[i]->storageFloat();

    if ( points[i]->energy() == 0.0 ) {
      if ( points[i]->use_sf() &&
           (points[i]->domain()->numNeighbors() == 6) ) {
        s4WSF.push_back(i);
      } else {
        s4.push_back(i);
      }
    } else if ( !((bias < split_val) || (bias > delete_val)) ) {
      continue;   

    // otherwise, should we split or delete this Point?
    } else if ( bias > 0.0 ) { // DELETE
      delP.push_back(i);
      ++num_deleted;
    } else { // SPLIT
      split.push_back(i);
      ++num_split;
    }
  }
  _optimized = false;

  //splitIntoFourWSF
  vector<int>::iterator iter = s4WSF.begin();
  while (iter != s4WSF.end()) {
    int i = *iter++;
    splitIntoFourWSF(i, points[i]);
  }

  //splitIntoFour
  iter = s4.begin();
  while (iter != s4.end()) {
    int i = *iter++;
    splitIntoFour(i, points[i]);
  }

  //splitAPoint
  iter = split.begin();
  while (iter != split.end()) {
    int i = *iter++;
    splitAPoint(i, points[i]);
  }

  //deleteAPoint
  sort(delP.begin(), delP.end());
  vector<int>::reverse_iterator riter = delP.rbegin();
  while (riter != delP.rend()) {
    int i = *riter++;
    deleteAPoint(i, points[i]);
  }

  cout << endl;
  cout << "Num split = " << num_split << "    Num delete = " << 
    num_deleted << "  " << points.size() << endl << endl;
  double seconds =  (double)(clock() - st)/(double)CLOCKS_PER_SEC;
  cout << "GlobalSurfaceEnergyNA::biasedOptimize iteration: "
       << seconds << " seconds." << endl;
}

//------------------------------------------------------------------------
// Function    : randomOptimize
// Description : optimze the system by randomly splitting and deleting
//               Points
//------------------------------------------------------------------------
void GlobalSurfaceEnergyNA::randomOptimize( svector<DynamicSurfacePoint*> 
                                            &points )
{
  _optimized = true;
  if ( points.empty() )
    return;

  _optimized = false;
  
  float ideal_local_energy = points[0]->domain()->idealEnergy();

  // sum up the global energy and ideal energy
  float energy=0.0, ideal_energy=0.0, tmp_energy;
  for ( int i = 0; i < points.size(); i++ )
  {  
    points[i]->computeEnergy( tmp_energy );
    energy += tmp_energy;
    ideal_energy += ideal_local_energy;
  }

  // check if the energy is zero
  if ( !energy )
  {
    splitPoints( points.size(), points );
    return;
  }

  // now, see if we need to split or kill any Points
  //
  // check if we are near this energy
  // --> adding one Point really works like increasing the 
  //     energy by the number of neighbors ... (i think?)
#ifdef THREE_D
  float energy_diff = (energy - ideal_energy) / 
                      (ideal_energy);
#endif

#ifdef TWO_D
  float energy_diff = (energy - ideal_energy) / 
                      (2.0*ideal_energy);
#endif

  if ( fabs(energy_diff) < 0.025 )
  {
    _optimized = true;
    cout << "Global Surf Energy OPTIMIZED" << endl;
    return;
  }
  
  // otherwise, we need to do something to drive energy to the ideal
  //   value

  if ( energy_diff > 0 )
  {
    // DELETE

    // first, determine approx how many particles we need to delete
    int num_delete =
      (int)(energy_diff * (float)points.size());

    // make sure we stay above the min number needed 
    //    
    // TODO : for now i have hardcoded the "minimum num points" to
    //        be 10!
    //
    num_delete = min( num_delete, (int)(points.size()-10) );

    if ( num_delete )
    {
      cout << "   Deleting " << num_delete << " particles...\n";
      deletePoints( num_delete, points );
    }
    else
      _optimized = true;
  }

  else 
  {
    // SPLIT
    
    // first, determine approx how many particles we need to add
    int num_add =
      (int)(-energy_diff * (float)points.size());

    if ( num_add )
    {
      cout << "   Splitting " << num_add << " particles...\n";
      splitPoints( num_add, points );
    }
    else
      _optimized = true;
  }
}


