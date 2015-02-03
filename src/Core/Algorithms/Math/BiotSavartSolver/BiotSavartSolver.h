
#ifndef CORE_ALGORITHMS_MATH_BIOTSAVARTSOLVER_H
#define CORE_ALGORITHMS_MATH_BIOTSAVARTSOLVER_H 1

//! Datatypes that the algorithm uses
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>

//! Base class for algorithm
#include <Core/Algorithms/Util/AlgoBase.h>

//! for Windows support
#include <Core/Algorithms/Fields/share.h>

namespace SCIRunAlgo {

using namespace SCIRun;

class SCISHARE BiotSavartSolverAlgo : public AlgoBase
{
  public:
    BiotSavartSolverAlgo()
    {
      // Number of processors to use
      add_int("num_processors",-1);
      istep = 0.0d;
    }
    
    //! Solve magnetic field via Biot-Savart piece-wise linear intergation
    bool run(FieldHandle& mesh, FieldHandle& coil, int outtype, MatrixHandle& outdata);

    //! For testing purpose, explicitly set the integration step
    void SetIntegrationStep(double step)
    {
    	istep = step;
    }

    //! For testing purpose, get the explicit integration step
    //! value less than 0 denotes no explicit step, auto-adjust step in effect 
    double GetIntegrationStep(void)
    {
    	return istep;
    }

private:
	double istep;

};

} // end namespace BiotSavartSolverAlgo

#endif
