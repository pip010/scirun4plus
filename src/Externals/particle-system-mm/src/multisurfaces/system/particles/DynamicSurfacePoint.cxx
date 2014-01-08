#include <cstdlib>
#include <iostream>
#include <system/domain/Domain.h>
#include <system/defines.h>
#include <system/particles/DynamicSurfacePoint.h>

#ifdef _WIN32
#pragma warning( disable : 4244 4305 4267 )
#endif

using namespace particle_sys;
using namespace std;

//#define DO_NORMAL_CHECK

//------------------------------------------------------------------------
// Function    : Constructors
// Description : 
//------------------------------------------------------------------------
DynamicSurfacePoint::DynamicSurfacePoint( Domain *d, ParticleSystem *sys, 
                                          const vector_type &pos)
  : SurfacePoint(d,sys,pos) 
{ 
  _angular_density = 0.0;
  _planar_seperation = 0.5;

  if ( d->insideDomain( pos ) )
    _moved_outside = false;
  else
    _moved_outside = true;

  _lambda = INITIAL_LAMBDA;
  _energy = 0.0;
  _force = _neighborhood_force = 0.0;

  _temporary = false;
  _using_sf = false;

  _storage_float = 0.0;
  _storage_int = 0;

  _max_stretch = MAX_VALUE;
}

DynamicSurfacePoint::DynamicSurfacePoint( const DynamicSurfacePoint &that )
  : SurfacePoint() 
{ *this = that; }

DynamicSurfacePoint::DynamicSurfacePoint() : SurfacePoint() 
{ 
  _angular_density = 0.0;
  _planar_seperation = 0.5;

   _moved_outside = false;

  _lambda = INITIAL_LAMBDA;
  _energy = 0.0;
  _force = _neighborhood_force = 0.0;

  _domain = NULL;
  _system = NULL;

  _temporary = false;

  _storage_float = 0.0;
  _storage_int = 0;

  _max_stretch = MAX_VALUE;
}

//------------------------------------------------------------------------
// Function    : computeEnergyAndForce
// Description : this is the big workhorse of this class! here is where
//               we go through the neighboring Points and compute the
//               local energy and force values
//
//               NOTE : this assumes a global radius of influence! use
//                      this function ONLY for doing the Gauss-Seidel
//                      updating scheme with a global radius
//------------------------------------------------------------------------
void DynamicSurfacePoint::computeEnergyForceYank( float &energy, 
                                                  vector_type &force, 
                                                  matrix_type &yank,
                                                  bool compute_energy,
                                                  bool compute_force,
                                                  bool compute_yank )
{
  const DynamicSurfacePoint *neighbor; // the neighboring Point
  float smoothing;  // smoothing term, based on normals

  // initialize the energy and force values
  clearEnergyForceYankBuffers();

  // initialize the Neighborhood and get the first neighbor
  _domain->determineNeighborhood( this );
  neighbor = _domain->nextNeighbor( this );

  // iterate through all of the neighboring points
  vector_type d_ij;
  while( neighbor )
  {
    // first, we want to eliminate any Points that are more than
    //   90 degrees away
#ifdef DO_NORMAL_CHECK
    if ( DotProduct(_normal, neighbor->normal()) <= 0.0 )
    {
      neighbor = _domain->nextNeighbor( this );
      continue;
    }
#endif

    // get the smoothing term of these two points, based upon their
    //   normals --> this reduces the effects of Points that are
    //   almost orthgonal to each other
    smoothing = smoothingCoefficient( _normal, neighbor->normal() );

    // compute the scaled vector between these two points

    // Does this distance use physical coordinates?  RTW - July 8, 2009
    distance( neighbor, d_ij );

    // initialize the energy
    _domain->initializeEnergy( d_ij );

    //
    // find the force and energy between the two Points
    //
    if ( compute_energy )
      _energy += smoothing*_domain->solveEnergy();
                                                 
    if ( compute_force )
      _force += smoothing*_domain->solveForce();

    if ( compute_yank )
      _yank += smoothing*_domain->solveYank();

    neighbor = _domain->nextNeighbor( this );
  }

  // check if we need to do a projection of the force
  if ( domain()->surface()->projectionVector() )
  {
    vector_type proj; 
    domain()->surface()->getProjectionVector( _position, proj );
    _force = DotProduct(_force, proj) * proj;
  }

  // return the energy and force buffer values
  energy = _energy;
  force = _force;
  yank = _yank;
}

void DynamicSurfacePoint::computeEnergy( float &energy )
{
  vector_type tmp_v;
  matrix_type tmp_m;
  computeEnergyForceYank( energy, tmp_v, tmp_m, true, false, false );
}

void DynamicSurfacePoint::computeForce( vector_type &force )
{
  float tmp_f;
  matrix_type tmp_m;
  computeEnergyForceYank( tmp_f, force, tmp_m, false, true, false );
}

void DynamicSurfacePoint::computeEnergyForce( float &energy, 
                                              vector_type &force )
{
  matrix_type tmp_m;
  computeEnergyForceYank( energy, force, tmp_m, true, true, false );
}

void DynamicSurfacePoint::computeYank( matrix_type &yank )
{
  float tmp_f;
  vector_type tmp_v;
  computeEnergyForceYank( tmp_f, tmp_v, yank, false, false, true );
}

//------------------------------------------------------------------------
// Function    : distance
// Description : compute the scaled vector between this particle and
//               one of its neighbors
//------------------------------------------------------------------------
void DynamicSurfacePoint::distance( const DynamicSurfacePoint *neighbor,
                                    vector_type &d_ij ) const
{
  // if using a sizing field
  if ( _using_sf )
  {
    // compute the average curvature at these two points
    //float ave_phi = 0.5*(sf() + neighbor->sf());
    float ave_phi = min(sf(), neighbor->sf());

    // compute the scaling factor
    float beta = _domain->normalizedIdealDistance() / (ave_phi+EPSILON);

    d_ij = beta * (_position - neighbor->position());

    return;
  }

  // compute the average curvature at these two points
  float ave_kappa = 0.5*(curvature_mag() + neighbor->curvature_mag());

  // compute the scaling factor
  float beta = stretch( ave_kappa ) / 
    (_planar_seperation/_domain->normalizedIdealDistance());

  //cout << beta << endl;

  d_ij = beta * (_position - neighbor->position());


  

  //matrix_type ave_kappa;
  //ave_kappa = 0.5*(curvature() + neighbor->curvature());

  //matrix_type beta;
  //stretch( ave_kappa, beta );
  //beta /= (_planar_seperation/NORMALIZED_IDEAL_DISTANCE);

  //d_ij = beta * (_position - neighbor->position());

}

//------------------------------------------------------------------------
// Function    : epsilonSample
// Description : 
//------------------------------------------------------------------------
void DynamicSurfacePoint::epsilonSample( float &closest_pt, float &sf )
{
  closest_pt = MAX_VALUE;
  sf = _surf_params._sf;

  const DynamicSurfacePoint *neighbor; // the neighboring Point

  // initialize the Neighborhood and get the first neighbor
  _domain->determineNeighborhood( this );
  neighbor = _domain->nextNeighbor( this );

  // iterate through all of the neighboring points
  vector_type d_ij;
  while( neighbor )
  {
    // first, we want to eliminate any Points that are more than
    //   90 degrees away
#ifdef DO_NORMAL_CHECK
    if ( DotProduct(_normal, neighbor->normal()) <= 0.0 )
    {
      neighbor = _domain->nextNeighbor( this );
      continue;
    }
#endif

    // compute the scaled vector between these two points
    d_ij = _position - neighbor->position();

    closest_pt = min( d_ij.length(), closest_pt );

    neighbor = _domain->nextNeighbor( this );
  }
}

//------------------------------------------------------------------------
// Function    : stretch
// Description : 
//------------------------------------------------------------------------
float DynamicSurfacePoint::stretch( float curvature ) const
{
  float s = 1.0+(angular_density()/(2.0*PI))*curvature*_planar_seperation;
  return min( s, _max_stretch);
}

void DynamicSurfacePoint::stretch( const matrix_type &curvature,
                                   matrix_type &s ) const
{
  s.identity();
  s += (angular_density()/(2.0*PI))*_planar_seperation*curvature;
}

//------------------------------------------------------------------------
// Function    : computeNeighborhoodAverageForce
// Description : compute the average force of the neighborhood -- the bool 
//               variable indicates whether we want to use the averaged 
//               force stored at each neighboring point (only useful after 
//               this function has been called once!)
//
//               TODO: should we add in the contribution from this point
//                     too?
//------------------------------------------------------------------------
void DynamicSurfacePoint::computeNeighborhoodAverageForce( 
  bool use_previous_averages )
{
  DynamicSurfacePoint *n;
  svector<DynamicSurfacePoint*> neighbors; // the neighboring Point
  
  // initialize the Neighborhood and get the first neighbor
  _domain->determineNeighborhood( this );
  n = _domain->nextNeighbor( this );

  // get the neighbors (need to get neighbors first because the
  //    computeForce() functions will reinitialize the domain!)
  while( n )
  {
    neighbors.push_back( n );
    n = _domain->nextNeighbor( this );
  }

  // average the force of the neighborhood
  vector_type f(0.0), f_tmp;
  for ( int i = 0; i < (int)neighbors.size(); i++ )
  {
    if ( !use_previous_averages )
      neighbors[i]->computeForce( f_tmp );
    else
      neighbors[i]->neighborhood_force();

    f += f_tmp;
  }

  // average the force
  if ( !neighbors.empty() )
    _neighborhood_force = f/(float)neighbors.size(); 
  else
    _neighborhood_force = 0.0;
}

//------------------------------------------------------------------------
// Function    : move
// Description : move the point to this new location, then do all the 
//               parameter  updating (the surface parameters and 
//               computation of the ideal energy at this new location)
//------------------------------------------------------------------------
void DynamicSurfacePoint::move( const vector_type &new_pos )
{
  // first, check if this new position is outside of the domain
  //   --> if so, keep old position and set flag!
  if ( !_domain->insideDomain( new_pos ) )
  {
    _moved_outside = true;
    return;
  }

  _moved_outside = false;


  // notify the Domain that the Point will be moving to a new 
  //   location (if this isn't a temporary point!)
  if ( !_temporary )
    _domain->movingToNewPosition( this, new_pos );

  // set the new position
  _position = new_pos;

  // update the surface parameters
  if ( !updateSurfaceParameters() )
  {
    _moved_outside = true;
    return;
  }

  computeRadius();
}

//------------------------------------------------------------------------
// Function    : computeRadius
// Description : 
//------------------------------------------------------------------------
void DynamicSurfacePoint::computeRadius()
{
//  //
//  // do a hacky curvature estimation here!!
//  //
//
//  const DynamicSurfacePoint *neighbor; // the neighboring Point
//
//  // initialize the Neighborhood and get the first neighbor
//  _domain->determineNeighborhood( this );
//  neighbor = _domain->nextNeighbor( this );
//
//  CotanWeightingFunction weight;
//  float summed_w=0.0, w_i;
//
//  // iterate through all of the neighboring points
//  vector_type r_ij, n_ij;
//  float r_ij_norm, summed_curvature=0.0;
//  int num_neighbors = 0;
//  while( neighbor )
//  {
//    // first, we want to eliminate any Points that are more than
//    //   90 degrees away
//#ifdef DO_NORMAL_CHECK
//    if ( DotProduct(_normal, neighbor->normal()) <= 0.0 )
//    {
//      neighbor = _domain->nextNeighbor( this );
//      continue;
//    }
//#endif
//
//    ++num_neighbors;
//    r_ij = neighbor->position() - _position;
//    r_ij_norm = r_ij.length();
//
//    n_ij = neighbor->normal() - _normal;
//
//    w_i = weight.weight( _influence_radius, r_ij_norm );
//    //w_i = 1.0;
//    summed_w += w_i; 
//
//    summed_curvature += 
//      w_i*pow(DotProduct( r_ij, n_ij ) / (r_ij_norm*r_ij_norm+EPSILON), 2); 
//    /*summed_curvature += 
//      w_i*n_ij.length() / (r_ij_norm*r_ij_norm+EPSILON);*/
//
//    neighbor = _domain->nextNeighbor( this );
//  }
//
//  _surf_params._curvature_mag = 
//    summed_curvature / (summed_w+EPSILON);
//  _surf_params._curvature_mag = sqrt( _surf_params._curvature_mag );


  // and set the radius of the Point to be half of the ideal seperation
  //   of this Point (because the ideal seperation is between Point
  //   centers!)
  
  // if using a sizing field
  if ( _using_sf )
  {
    _radius  = 0.5*sf();
    _radius2 = 0.5*sf();
    return;
  }

  _radius  = 0.5*_planar_seperation / stretch(curvature_mag());
  _radius2 = 0.5*_planar_seperation / stretch(curvature_mag());
  

  //_radius = 0.5*_planar_seperation / stretch(kappa1());
  //_radius2 = 0.5*_planar_seperation / stretch(kappa2());
}

//------------------------------------------------------------------------
// Function    : offset
// Description : offset this point along it's normal direction
//------------------------------------------------------------------------
void DynamicSurfacePoint::offset( float offset_amt )
{
  vector_type new_position = _position + offset_amt*_normal;

  // first, check if this new position is outside of the domain
  //   --> if so, keep old position and set flag!
  if ( !_domain->insideDomain( new_position ) )
  {
    _moved_outside = true;
    return;
  }
  _moved_outside = false;

  // notify the Domain that the Point will be moving to a new 
  //   location 
  _domain->movingToNewPosition( this, new_position );

  // set the new position
  _position = new_position;
}

//------------------------------------------------------------------------
// Function    : smoothingCoefficient
// Description : returns a smoothing value based on the angle between
//               two Points' neighbors
//------------------------------------------------------------------------
float DynamicSurfacePoint::smoothingCoefficient( const vector_type& p_i_normal,
                                                 const vector_type& p_j_normal )  
{
#ifndef DO_NORMAL_CHECK
  return 1.0;
#endif

  // get the cosine of the angle between the two particles' normals
  float cosine = DotProduct(p_i_normal,p_j_normal)/
      (p_i_normal.length()*p_j_normal.length() + EPSILON);

  // define the float region cut-off in terms of the cosine
  float flat_cutoff = 0.156; // approx 90% of the region

  // the flat region
  if ( cosine >= flat_cutoff ) return 1.0;

  // the feathered region
  return ( cos((flat_cutoff-cosine)/flat_cutoff* (PI/2.0)) ); 
}

//------------------------------------------------------------------------
// Function    : operator =
// Description : the assignment operator
//------------------------------------------------------------------------
const DynamicSurfacePoint& DynamicSurfacePoint::operator = (
  const DynamicSurfacePoint &that)
{
  // copy the SurfacePoint variables
  (*this).surfacePointAssignmentOp( that );

  _angular_density = that._angular_density;
  _planar_seperation = that._planar_seperation;
  _moved_outside = that._moved_outside;
  _lambda = that._lambda;
  _max_stretch = that._max_stretch;

  _storage_float = 0.0;
  _storage_int = 0;

  _using_sf = that._using_sf;

  return *this;
}

void DynamicSurfacePoint::dynamicSurfacePointAssignmentOp(
  const DynamicSurfacePoint &that) 
{ 
  *this = that; 
}

//------------------------------------------------------------------------
// Function    : operator << 
// Description : output function for the DynamicSurfacePoint class
//------------------------------------------------------------------------
ostream& operator << ( ostream& outs, const DynamicSurfacePoint& source ) 
{
  outs << (SurfacePoint)source << endl;
  
  return outs; 
}



