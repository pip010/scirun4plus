
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
    {}
    
    //! Convert data into a matrix
    bool run(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& outdata);

  protected:
    //! Biot-Savart Contour Piece-wise integration
    bool IntegrateBiotSavart(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& outdata);
};

} // end namespace BiotSavartSolverAlgo

#endif
