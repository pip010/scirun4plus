
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
      istep = 0.0;
      tfactor = 0u;
    }
    
    //! Solve magnetic field via Biot-Savart piece-wise linear intergation
    bool run(FieldHandle& mesh, FieldHandle& coil, int outtype, MatrixHandle& outdata);

    //! For testing purpose, explicitly set the integration step
    //! thus auto adaptation steps will be skiped (good for testing/experiments)
    void SetIntegrationStep(double step)
    {
    	istep = step;
    }

    //! For testing purpose, get the explicit integration step
    //! value of 0 denotes no explicit step that is auto adaptation is active
    double GetIntegrationStep(void)
    {
    	return istep;
    }

    //! For testing purposes, explicitly set the number of threads to multithread
    //! value of 1 or 0 will disable multithread
    //! value of 4 is recommended for typical quadcore CPUs
    void SetThreadingFactor(unsigned int factor)
    {
    	tfactor = factor;
    }

    //! 
    unsigned int GetThreadingFactor(void)
    {
    	return tfactor;
    }

private:
	double istep;
	unsigned int tfactor;
};

} // end namespace BiotSavartSolverAlgo

#endif
