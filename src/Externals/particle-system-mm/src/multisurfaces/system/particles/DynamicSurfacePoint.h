//----------------------------------------------------------------------//
// FILE        : DynamicSurfacePoint.h                                                
// DESCRIPTION : The DynamicSurfacePoint class extends the SurfacePoint 
//               and Point class, and provides additional dynamic 
//               characteristics to Points associated with Surfaces. The
//               Points now have a notion of energy and force due to
//               their neighboring Points. DynamicSurfacePoints can 
//               also be moved. 
//----------------------------------------------------------------------//

#ifndef __DYNAMIC_SURFACE_POINT_H__ 
#define __DYNAMIC_SURFACE_POINT_H__ 

#include <system/systemExports.h>

#include <features/mtxlib.h>
#include <features/svector.h>
#include <system/defines.h>
#include <system/particles/SurfacePoint.h>

namespace particle_sys 
{ 
  template <class T>
  class Neighborhood;

  class System_SHARE DynamicSurfacePoint : public SurfacePoint
  {
  public:
    DynamicSurfacePoint(Domain *d, ParticleSystem *sys, const vector_type &pos);
    DynamicSurfacePoint(const DynamicSurfacePoint &that);
    DynamicSurfacePoint();
    virtual ~DynamicSurfacePoint() {}

    const DynamicSurfacePoint& operator = (const DynamicSurfacePoint &that);

    //--------------------------------------------------------------------
    // energy and force calculations --> two functions are provided so 
    //   that the energy or force only can be computed 
    //
    //   someone on the outside should never have to call 
    //      computeEnergyAndForce() with the bool variables!
    virtual void computeEnergyForceYank(float &energy, vector_type &force,
                                        matrix_type &yank,
                                        bool compute_energy=true, 
                                        bool compute_force=true,
                                        bool compute_yank=true);
    void computeEnergy(float &energy);
    void computeForce(vector_type &force);
    void computeEnergyForce(float &energy, vector_type &force);
    void computeYank(matrix_type &yank);

    inline void clearEnergyForceYankBuffers() 
    { _energy = 0.0; _force = 0.0; _yank = 0.0;}

    //--------------------------------------------------------------------
    // compute the scaled vector between this particle and one of its
    //    neighbors
    void distance(const DynamicSurfacePoint *neighbor, vector_type &d_ij) 
      const;

    //--------------------------------------------------------------------
    // compute the average force of the neighborhood -- the bool variable
    //   indicates whether we want to use the averaged force stored at
    //   each neighboring point (only useful after this function has been
    //   called once!)
    void computeNeighborhoodAverageForce(bool use_previous_averages=false);

    //--------------------------------------------------------------------
    // move the point to this new location, then do all the parameter 
    //   updating (the surface parameters and computation of the ideal
    //   energy at this new location)
    void move(const vector_type &new_pos);

    void offset(float offset_amt);

    //--------------------------------------------------------------------
    // resets the Point's lambda to the initial value (this value is set
    //   in defines.h)
    inline void resetLambda() { _lambda = INITIAL_LAMBDA; };

    //--------------------------------------------------------------------
    // setter and getter functions to communticate the System's angular 
    //   density and planar seperation values
    inline float angular_density() const { return _angular_density; };
    inline float planar_seperation() const { return _planar_seperation; };
    inline bool  moved_outside() const { return _moved_outside; };
    inline float lambda() const { return _lambda; };
    inline float energy() const { return _energy; };
    inline vector_type force() const { return _force; };
    inline matrix_type yank() const { return _yank; };

    inline void angular_density(float ad) 
    { _angular_density = ad; resetLambda(); };
    inline void planar_seperation(float ps) 
    { _planar_seperation = ps; resetLambda(); };
    inline void lambda(float l) { _lambda = l; };

    inline const vector_type& neighborhood_force() const
    { return _neighborhood_force; };

    inline custom_class::svector<DynamicSurfacePoint*>* neighbors()
    { return &_neighbors; };

    //--------------------------------------------------------------------
    // access to the _temporary variable, which is used for creating new
    //   DynamicSurfacePoints for potential addition to the System
    inline void temporary(bool t) { _temporary = t; };
    inline bool temporary() { return _temporary; };

    //--------------------------------------------------------------------
    // storage variables that can be used 
    inline void storageFloat(float sf) 
    { _storage_float = sf; };
    inline float storageFloat() const
    { return _storage_float; };
    inline void storageInt(int si) 
    { _storage_int = si; };
    inline int storageInt() const
    { return _storage_int; };

    //--------------------------------------------------------------------
    // compute the influence radius of the particle
    virtual void computeRadius();

    //--------------------------------------------------------------------
    // compute the minimum seperation based on a percentage of the
    //   planar seperation (the max allowed!)
    inline void minSeperationPercentage(float percentage)
    { _max_stretch = 1.0/(percentage+EPSILON); };


    inline float radius2() const { return _radius2; };

    inline void using_sf() { _using_sf = true; };
    inline bool use_sf() const { return _using_sf; };

    void epsilonSample(float &closest_pt, float &sf);

  protected:
    //--------------------------------------------------------------------
    // this will return the stretch of space related to this curvature
    float stretch(float curvature) const;
    void stretch(const matrix_type &curvature, matrix_type &s) const;

    void dynamicSurfacePointAssignmentOp(const DynamicSurfacePoint &that); 

    float _planar_seperation, _max_stretch;
    bool  _moved_outside;

    float _lambda, // for adaptive dt and LM integration
          _energy;
    vector_type _force, _neighborhood_force;
    matrix_type _yank;

    bool _temporary;  // is this a temp point? (for splitting)

    float _storage_float;
    int _storage_int;

    float _radius2;

    bool _using_sf;

    custom_class::svector<DynamicSurfacePoint*> _neighbors;

    //--------------------------------------------------------------------
    // this is a smoothing function based on the angle between to
    //   Points' normals
    float smoothingCoefficient(const vector_type& p_i_normal, 
                               const vector_type& p_j_normal);

  };  

} // namespace particle_sys 

// print function
std::ostream& operator << (std::ostream& outs, 
                           const particle_sys::DynamicSurfacePoint& source);
#endif // __DYNAMIC_SURFACE_POINT_H__ 

