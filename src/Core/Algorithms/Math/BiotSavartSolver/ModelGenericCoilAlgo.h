
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
  	std::vector<Vector> ConcatPointsForCurve(std::vector<Vector>& points1, std::vector<Vector>& points2);
  	std::vector<size_t> ConcatIndicesForCurve(std::vector<size_t>& indices1,std::vector<size_t>& indices2);
  	void GenerateCircleContour(std::vector<Vector>& points, std::vector<size_t>& indices, Vector pos, double r, size_t nsegments);
  	void GenerateFigure8ShapedCoil(std::vector<Vector>& points, std::vector<size_t>& indices, double r, double d, size_t nsegments);
};

} // end namespace BiotSavartSolverAlgo

#endif