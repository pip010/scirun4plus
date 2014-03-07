#ifndef __SCALAR_FIELD_H__
#define __SCALAR_FIELD_H__

#include "mtxlib.h"
#include "ScalarFieldParameters.h"
#include <iostream>


#define EPSILON (1.0e-12)


/************************************************************************/
//                   The Scalar Field Base Class                        //
/************************************************************************/
template <int dim> class ScalarField
{
public:
  ScalarField() {};
  virtual ~ScalarField(){};

  //---------------------------------------------//
  //           virtual functions                 //
  //---------------------------------------------//
     
  //--------------------------------------------------------------------
  // functions for updating scalar field values at particle locations
  virtual void computeScalarFieldParams(const vec<dim> &pos,
                                        ScalarFieldParams<dim> &params,
                                        bool computeHessian=false) 
    const = 0;

  matrix<dim,dim> hessian(const vec<dim> &v)
  { ScalarFieldParams<dim> p;
    computeScalarFieldParams( v, p, true );
    return p._Fxx; };
  vec<dim> gradient(const vec<dim> &v)
    { ScalarFieldParams<dim> p;
    computeScalarFieldParams( v, p, false );
    return p._Fx; };  
  mtx_datatype value(const vec<dim> &v)
  { ScalarFieldParams<dim> p;
    computeScalarFieldParams( v, p, false );
    return p._F; }; 

  virtual bool checkBounds(const vec<dim> &pos) const {return(true);}
  virtual vec<dim> boundingBoxLow() const
  { return(vec<dim>(0.0f)); };
  virtual vec<dim> boundingBoxHigh() const
  { std::cout << "got incorrect call to bbh" << std::endl; return(vec<dim>(0.0f)); }

  virtual int xdim () const {return(0);}
  virtual int ydim () const {return(0);}
  virtual int zdim () const {return(0);}

  //---------------------------------------------//
  //           surface functions                 //
  //---------------------------------------------//

  // compute the curvature at a particle based on gordon's
  //   vis03 paper
  mtx_datatype computeCurvature(ScalarFieldParams<dim> &params) const;

protected:

  mtx_datatype computeCurvatureMagnitude(const matrix<dim,dim> &G,
                                  ScalarFieldParams<dim> &params) const;

};

//------------------------------------------------------------------------
// Function    : computeCurvature()
// Description : compute the curvature at the point pos -- uses the algo
//               from gordon's vis03 paper
//------------------------------------------------------------------------
template <int dim>
mtx_datatype ScalarField<dim>::computeCurvature( ScalarFieldParams<dim> &params )
  const
{
  // get the normalized gradient at this point on the surface
  // TODO: i'm not sure why this is negative -- will have to look into 
  //       this later!!!!
  vec<dim> n = -params._Fx;
  n.normalize();

  // get the tangent plane projection matrix
  matrix<dim,dim> P; P.identity();
  P -= DirectProduct( n, n ); 

  // store the curvature matrix that has been projected into the 
  //   tangent plane
  matrix<dim,dim> G = (P*params._Fxx) / 
    (params._Fx.length()+EPSILON); // positive bc inside is lower values!

  // isolate the curvature values
  G *= P; 
  params._curvature = G;
  params._curvature_mag = computeCurvatureMagnitude( G, params );

  return params._curvature_mag;
}



#endif // __SCALAR_FIELD_H__
