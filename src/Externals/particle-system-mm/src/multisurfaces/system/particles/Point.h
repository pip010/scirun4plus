//----------------------------------------------------------------------//
// FILE        : Point.h                                                
// DESCRIPTION : The Point class is a basic object that stores a        
//               minimal amount of info of a point in space. This is the
//               core of the dynamic particles, and is used for point
//               based rendering.
//----------------------------------------------------------------------//

#ifndef __POINT_H__ 
#define __POINT_H__ 

#include <system/systemExports.h>

#include <features/mtxlib.h>
#include <system/defines.h>

#ifdef _WIN32
#pragma warning( disable : 4305 )
#endif

namespace particle_sys 
{ 
  class Domain;
  class ParticleSystem;

  class System_SHARE Point
  {
  public:
    Point(Domain *d, ParticleSystem *sys) 
    { _domain=d; _system=sys; _radius=1.0; };
    Point(Domain *d, ParticleSystem *sys, const vector_type &pos): 
      _position(pos), _normal(pos)                     
    { _domain=d; _system=sys; _radius=0.1; _normal.normalize(); };
    Point(const vector_type &pos, const vector_type &norm, float r)
    { _position = pos; _normal = norm; _radius = r; };
    Point(const Point &that) { *this = that; };
    ~Point() {};

    const Point& operator = (const Point &that)
    {
      _position = that._position;
      _normal = that._normal;
      _radius = that._radius;
      _domain = that._domain;
      _system = that._system;

      return *this;
    };

    //--------------------------------------------------------------------
    // getter functions for the basic Point variables
    inline const vector_type& position() const { return _position; };
    inline const vector_type& normal() const { return _normal; };
    inline float radius() const { return _radius; };
    inline Domain* domain() const { return _domain; };
    inline ParticleSystem* system() const { return _system; };

    inline void radius(float r) { _radius = r; };
    inline void normal(vector_type &n) { _normal = n; };
    inline void domain(Domain *d) { _domain = d; };

  protected:
    vector_type _position, _normal;
    float   _radius;
    Domain *_domain;
    ParticleSystem *_system;

    Point() { _domain=NULL; _system=NULL; _radius=0.0; };
    void pointAssignmentOp(const Point &that) { *this = that; };
  };  

} // namespace particle_sys 

// print function
std::ostream& operator << (std::ostream& outs, 
                           const particle_sys::Point& source);
#endif // __POINT_H__ 
