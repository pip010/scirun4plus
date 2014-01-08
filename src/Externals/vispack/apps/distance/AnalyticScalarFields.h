#ifndef __ANALYTIC_SCALAR_FIELDS_H__
#define __ANALYTIC_SCALAR_FIELDS_H__

#include "ScalarField.h"
#include "mtxlib.h"
#include <iostream>

/************************************************************************/
//   These analytic fields are always centered at the origin.           //
/************************************************************************/

/************************************************************************/
//                            SPHERE  - center of 1,1,1, box                                   //
/************************************************************************/
class SpherePlane : public ScalarField<3>
{
public:
  SpherePlane(mtx_datatype radius, mtx_datatype z) { _radius = radius; _z = z;}
  ~SpherePlane(){};

  void computeScalarFieldParams(const vec<3> &pos, 
                                ScalarFieldParams<3> &params,
                                bool computeHessian=false) const;

  virtual vec<3> boundingBoxLow() const
  { return(vec<3>(0.0f)); };
  virtual vec<3> boundingBoxHigh() const
  { std::cout << "got correct call to bbh" << std::endl; return(vec<3>(1.0f, 1.0f, 1.0f)); }

private:
  mtx_datatype _radius;
  mtx_datatype _z;
};



/************************************************************************/
//                            SPHERE                                    //
/************************************************************************/
template <int dim> class Sphere : public ScalarField<dim>
{
public:
  Sphere(mtx_datatype radius) { _radius = radius; };
  ~Sphere(){};

  void computeScalarFieldParams(const vec<dim> &pos, 
                                ScalarFieldParams<dim> &params,
                                bool computeHessian=false) const;

private:
  mtx_datatype _radius;
};


/************************************************************************/
//                            TORUS                                     //
/************************************************************************/
template <int dim> class Torus : public ScalarField<dim>
{
public:
  Torus(mtx_datatype r, mtx_datatype R) { _r = r; _R = R; }; 
  ~Torus(){};

  void computeScalarFieldParams(const vec<dim> &pos, 
                                ScalarFieldParams<dim> &params,
                                bool computeHessian=false) const;
      
private:
  mtx_datatype _r, _R;
};


/************************************************************************/
//                        QUARTIC CUBE                                  //
/************************************************************************/
template <int dim> class QuarticCube : public ScalarField<dim>
{
public:
  QuarticCube(mtx_datatype radius) { _radius = radius; };
  ~QuarticCube(){};

  void computeScalarFieldParams(const vec<dim> &pos, 
                                ScalarFieldParams<dim> &params,
                                bool computeHessian=false) const;

private:
  mtx_datatype _radius;
};


#endif // __ANALYTIC_SCALAR_FIELDS_H__
