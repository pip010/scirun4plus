#include <cstdlib>
#include <system/particles/DynamicSurfacePoint.h>
#include <system/domain/Intersection.h>

using namespace std;
using namespace particle_sys;

//------------------------------------------------------------------------
// Function    : constructor and destructor
// Description : 
//------------------------------------------------------------------------
Intersection::Intersection( IOScalarField *io, 
                            int num_indicators, 
                            int *indicators ) :Surface()
{
  // get the list of indicator function indices
  _num_indicators = num_indicators;
  _indicators = new int[_num_indicators];
  for ( int i = 0; i < _num_indicators; i++ )
    _indicators[i] = indicators[i];

  _io_field = io;
}

Intersection::~Intersection()
{
  delete [] _indicators;
}

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : 
//------------------------------------------------------------------------
bool Intersection::computeSurfacePointParams( const vector_type &pos, 
                                              SurfacePointParams 
                                              &params,
                                              bool computeHessian ) const
{  
  float E=0.0, sf=MAX_VALUE;
  vector_type Ex(0.0);

  for ( int i = 0; i < _num_indicators; i++ )
  {
    _io_field->setMainIndicatorFunction( _indicators[i] );
    _io_field->computeSurfacePointParams( pos, params, computeHessian );

    E  += (params._F*params._F);
    Ex += (2.0*params._F*params._Fx);
    sf = min( params._sf, sf );
  }

  params._F = E;
  params._Fx = Ex;
  params._sf = sf;

  return true;
}

//------------------------------------------------------------------------
// Function    : projectMotion()
// Description : 
//------------------------------------------------------------------------
void Intersection::projectMotion( const vector_type &pos,
                                  const vector_type &n,
                                  vector_type &m ) const
{
  SurfacePointParams params;
  vector_type n1, n2, n3, P;
  switch ( _num_indicators )
  {
  case 2:
    _io_field->setMainIndicatorFunction( _indicators[0] );
    _io_field->computeSurfacePointParams( pos, params, false );

    n1 = params._Fx;
    //n1.normalize();

    _io_field->setMainIndicatorFunction( _indicators[1] );
    _io_field->computeSurfacePointParams( pos, params, false );

    n2 = params._Fx;
    //n2.normalize();

    P = n2-n1;
    P.normalize();

    m -= DotProduct( P, m ) * P;
    break;

#ifdef THREE_D
  case 3:
    _io_field->setMainIndicatorFunction( _indicators[0] );
    _io_field->computeSurfacePointParams( pos, params, false );

    n1 = params._Fx;
    n1.normalize();

    _io_field->setMainIndicatorFunction( _indicators[1] );
    _io_field->computeSurfacePointParams( pos, params, false );

    n2 = params._Fx;
    n2.normalize();

    _io_field->setMainIndicatorFunction( _indicators[2] );
    _io_field->computeSurfacePointParams( pos, params, false );

    n3 = params._Fx;
    n3.normalize();

    P = 0.0;
    P += CrossProduct( n1, n2 );
    P += CrossProduct( n2, n3 );
    P += CrossProduct( n3, n1 );
    
    P.normalize();

    

    //cout << n1 << endl;
    //cout << n2 << endl;
    //cout << n2 << endl << endl;

    m = DotProduct( m, P ) * P;

    break;
#endif

  }
}

//------------------------------------------------------------------------
// Function    : computeSurfacePointParams()
// Description : 
//------------------------------------------------------------------------
bool Intersection::projectOntoSurface( DynamicSurfacePoint *point,
                                       vector_type &proj ) const
{
  return false;

  // compute the smoothing factor, based on how close Eij is to zero
  float E_threshold = 0.001;
  
  float s, E=point->F();
  if ( E >= E_threshold )
    s = 1.0;
  else
    s = 0.5*(cos( ((E_threshold-E) / E_threshold)*PI )+1.0);

  // add in some fraction of the E projection
  proj = s*(-E)/ ((point->Fx()).length()+EPSILON) * point->normal();

  // if near the junction, we need to regularize
  if ( s != 1.0 )
  {
    //cout << s << "  " << E << endl;

    SurfacePointParams params;
    vector_type pv;
    float f1, f2, f3;
    vector_type fx1, fx2, fx3;
    switch ( _num_indicators )
    {
    case 2:
      // in this situation, near the 2-junction grad.f1 and grad.f2
      //   will be approximately equal and opposite, so we can take a
      //   short-cut and just project with f1
      _io_field->setMainIndicatorFunction( _indicators[0] );
      _io_field->computeSurfacePointParams( point->position(), 
                                        params, false );

      pv = (-params._F)/(params._Fx).lengthSqr() * params._Fx;
      break;

    case 3:
      // in this situation, near the 3-junction the gradiants should
      //   all be seperated by an angle of ~120, so the strategy is 
      //   to first project onto the plane of 2 materials, and then
      //   project down onto the line of the 3 materials
      _io_field->setMainIndicatorFunction( _indicators[0] );
      _io_field->computeSurfacePointParams( point->position(), 
                                        params, false );
      f1 = params._F;
      fx1 = params._Fx;

      _io_field->setMainIndicatorFunction( _indicators[1] );
      _io_field->computeSurfacePointParams( point->position(), 
                                        params, false );
      f2 = params._F;
      fx2 = params._Fx;

      _io_field->setMainIndicatorFunction( _indicators[2] );
      _io_field->computeSurfacePointParams( point->position(), 
                                        params, false );
      f3 = params._F;
      fx3 = params._Fx;


      pv = -0.5*( (f1/fx1.lengthSqr())*fx1 - (f2/fx2.lengthSqr())*fx2 ) -
        (f3/fx3.lengthSqr())*fx3;
      break;

    case 4:
      cout << "4-junction material projection not implemented!!!" << endl;
      break;
    }

    // do an average between old way and new way
    proj += (1.0-s)*pv;
  }

  return true;
}

