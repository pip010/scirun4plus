
#ifndef CORE_ALGORITHMS_MATH_MODELGENERICCOIL_H
#define CORE_ALGORITHMS_MATH_MODELGENERICCOIL_H 1

//! Datatypes that the algorithm uses
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Geometry/Vector.h>

//! Base class for algorithm
#include <Core/Algorithms/Util/AlgoBase.h>

//#include <Core/Algorithms/Converter/MatrixToField.h>

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

    struct Args
    {
    	double wireCurrent;
    	double coilRadius;
    	double coilDistance;
    	double coilSegments;
    	int type;
    	
    	inline bool operator==(const Args& rhs)
    	{ 
			return wireCurrent == rhs.wireCurrent && 
			coilRadius == rhs.coilRadius && 
			coilDistance == rhs.coilDistance && 
			coilSegments == rhs.coilSegments && 
			type == rhs.type;
    	}
		inline bool operator!=( const Args& rhs){return !( (*this) == rhs );}
    };
    
    //! Convert data into a matrix
    bool run(FieldHandle& mesh, MatrixHandle& params, Args& args); 
};

} // end namespace BiotSavartSolverAlgo

#endif
