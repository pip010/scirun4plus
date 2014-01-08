#ifndef __ENERGY_H__
#define __ENERGY_H__

#include <system/systemExports.h>

#include <system/defines.h>
#include <features/mtxlib.h>

namespace particle_sys
{
  /**********************************************************************/
  // Energy is the kernel of the motion and life of the particules' 
  // System. By defining, the neighborhood and a energy's value 
  // depending on the distance between particles and the way to compute 
  // the neighbor
  /**********************************************************************/


  /**********************************************************************/
  //                     The Energy Base Class                          //
  /**********************************************************************/
  class System_SHARE Energy
  {
      
  public:
    // constructor
    Energy(float nid) {_normalized_ideal_distance = nid;};
    virtual ~Energy(){}

    //---------------------------------------------//
    //           virtual functions                 //
    //---------------------------------------------//

    virtual void initializeParameters(const vector_type &d_ij) = 0;

    // solve for the force and energy of the specific force function
    virtual float solveEnergy() = 0;
    virtual vector_type solveForce() = 0;
    virtual matrix_type solveYank() = 0;

    //---------------------------------------------//
    //           base clase functions              //  
    //---------------------------------------------//
    inline float idealEnergy() const
    { return _ideal_energy; };

  protected:
    float _ideal_energy;
    float _e, _f, _y; 

    float _normalized_ideal_distance;

  };

  /************************************************************************/
  //                          Cotan Energy                                //
  /************************************************************************/    
  class System_SHARE CotanEnergy : public Energy
  {
  public:
    CotanEnergy(float nid);
    ~CotanEnergy(){};

    void initializeParameters(const vector_type &d_ij);

    float solveEnergy();
    vector_type solveForce();
    matrix_type solveYank();

  private:
    float _d;
    vector_type _d_ij_unit; 
  };

  /************************************************************************/
  //                          Radial Energy                                //
  /************************************************************************/    
  class System_SHARE RadialEnergy : public Energy
  {
  public:
    RadialEnergy(float nid);
    ~RadialEnergy(){};

    void initializeParameters(const vector_type &d_ij);

    float solveEnergy();
    vector_type solveForce();
    matrix_type solveYank(){matrix_type tmp(0.0); return tmp;};
	

  private:
    float _d;
    vector_type _d_ij_unit; 
  };


}

#endif
