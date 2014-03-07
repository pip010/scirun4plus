#ifndef __KERNEL_H__
#define __KERNEL_H__

#include "mtxlib.h"

class Kernel
{
public:
  Kernel(){};
  virtual ~Kernel(){};

  virtual float w(const vec<4> &p, const vec<4> &tau)
  {std::cout << "Kernel::PROBLEM\n"; return 0.0;};
  virtual float w(const vec<4> &p, const vec<3> &tau)
  {std::cout << "Kernel::PROBLEM\n"; return 0.0;};
  virtual float w(const vec<4> &p, const vec<2> &tau)
  {std::cout << "Kernel::PROBLEM\n"; return 0.0;};

protected:

};

class CatmullRom : public Kernel
{
public:
  CatmullRom() : Kernel()
  { M.set( vec<4>(-1.0, 2.0,-1.0, 0.0),
           vec<4>( 3.0,-5.0, 0.0, 2.0),
           vec<4>(-3.0, 4.0, 1.0, 0.0),
           vec<4>( 1.0,-1.0, 0.0, 0.0) ); M /= 2.0; };
  ~CatmullRom() {};

  inline float w(const vec<4> &p, const vec<4> &tau)
  { return DotProduct(p*M,tau); };

private:
  matrix<4,4> M;

};

class CatmullRomD : public Kernel
{
public:
  CatmullRomD() : Kernel()
  { M.set( vec<3>(-3.0,  4.0,-1.0),
           vec<3>( 9.0,-10.0, 0.0),
           vec<3>(-9.0,  8.0, 1.0),
           vec<3>( 3.0, -2.0, 0.0) ); M /= 2.0; };
  ~CatmullRomD() {};

  inline float w(const vec<4> &p, const vec<3> &tau)
  { return DotProduct(p*M,tau); };

private:
  matrix<4,3> M;

};

class CatmullRomDD : public Kernel
{
public:
  CatmullRomDD() : Kernel()
  { M.set( vec<2>(-3.0, 2.0),
           vec<2>( 9.0,-5.0),
           vec<2>(-9.0, 4.0),
           vec<2>( 3.0,-1.0) ); };
  ~CatmullRomDD() {};

  inline float w(const vec<4> &p, const vec<2> &tau)
  { return DotProduct(p*M,tau); };

private:
  matrix<4,2> M;

};

class CubicBSpline : public Kernel
{
public:
  CubicBSpline() : Kernel()
  { M.set( vec<4>(-1.0, 3.0,-3.0, 1.0),
           vec<4>( 3.0,-6.0, 0.0, 4.0),
           vec<4>(-3.0, 3.0, 3.0, 1.0),
           vec<4>( 1.0, 0.0, 0.0, 0.0) ); M /= 6.0; };
  ~CubicBSpline() {};

  inline float w(const vec<4> &p, const vec<4> &tau)
  { return DotProduct(p*M,tau); };

private:
  matrix<4,4> M;

};

class CubicBSplineD : public Kernel
{
public:
  CubicBSplineD() : Kernel()
  { M.set( vec<3>(-1.0, 2.0,-1.0),
           vec<3>( 3.0,-4.0, 0.0),
           vec<3>(-3.0, 2.0, 1.0),
           vec<3>( 1.0, 0.0, 0.0) ); M /= 2.0; };
  ~CubicBSplineD() {};

  inline float w(const vec<4> &p, const vec<3> &tau)
  { return DotProduct(p*M,tau); };

private:
  matrix<4,3> M;

};

class CubicBSplineDD : public Kernel
{
public:
  CubicBSplineDD() : Kernel()
  { M.set( vec<2>(-1.0, 1.0),
           vec<2>( 3.0,-2.0),
           vec<2>(-3.0, 1.0),
           vec<2>( 1.0, 0.0) ); };
  ~CubicBSplineDD() {};

  inline float w(const vec<4> &p, const vec<2> &tau)
  { return DotProduct(p*M,tau); };

private:
  matrix<4,2> M;

};

#endif // __KERNEL_H__
