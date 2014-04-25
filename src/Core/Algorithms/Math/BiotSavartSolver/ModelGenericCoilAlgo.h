
#ifndef CORE_ALGORITHMS_MATH_MODELGENERICCOIL_H
#define CORE_ALGORITHMS_MATH_MODELGENERICCOIL_H 1

//! Datatypes that the algorithm uses
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Geometry/Vector.h>

//! Base class for algorithm
#include <Core/Algorithms/Util/AlgoBase.h>

#include <Core/Algorithms/Converter/MatrixToField.h>

//! for Windows support
#include <Core/Algorithms/Fields/share.h>

#include <vector>

namespace SCIRunAlgo {

using namespace SCIRun;

class SCISHARE ModelGenericCoilAlgo : public AlgoBase
{
  public:
    ModelGenericCoilAlgo()
    {}
    
    //! Convert data into a matrix
    bool run(FieldHandle& mesh, MatrixHandle& params); 

  private:
  	std::vector<Vector>	GenerateCircleContour(Vector pos,double r,uint nsegments);
};

} // end namespace BiotSavartSolverAlgo

#endif