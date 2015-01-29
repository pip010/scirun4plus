
#ifndef CORE_ALGORITHMS_TMS_MODELGENERICCOIL_H
#define CORE_ALGORITHMS_TMS_MODELGENERICCOIL_H 1

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

    class SCISHARE ModelTMSCoilSingleAlgo : public AlgoBase
    {
      public:
        ModelTMSCoilSingleAlgo()
        {}

        struct Args
        {
        	double wireCurrent;
        	double coilRadius;
            double coilDistanceOuter;
        	size_t coilLevelDetails;
        	int type;
        	
        	inline bool operator==(const Args& rhs)
        	{ 
    			return wireCurrent == rhs.wireCurrent && 
    			coilRadius == rhs.coilRadius && 
    			coilLevelDetails == rhs.coilLevelDetails &&
                coilDistanceOuter == rhs.coilDistanceOuter &&
    			type == rhs.type;
        	}
    		inline bool operator!=( const Args& rhs){return !( (*this) == rhs );}
        };
        
        //! Generate the coil geom
        bool run(FieldHandle& mesh, Args& args); 
    };


    class SCISHARE ModelTMSCoilSpiralAlgo : public AlgoBase
    {
      public:
        ModelTMSCoilSpiralAlgo()
        {}

        struct Args
        {
            double wireCurrent;
            size_t wireLoops;
            double coilRadiusInner;
            double coilRadiusOuter;
            double coilDistanceOuter;
            size_t coilLevelDetails;
            int type;
            
            inline bool operator==(const Args& rhs)
            { 
                return wireCurrent == rhs.wireCurrent && 
                wireLoops == rhs.wireLoops &&
                coilRadiusInner == rhs.coilRadiusInner && 
                coilRadiusOuter == rhs.coilRadiusOuter && 
                coilLevelDetails == rhs.coilLevelDetails &&
                coilDistanceOuter == rhs.coilDistanceOuter &&
                type == rhs.type;
            }
            inline bool operator!=( const Args& rhs){return !( (*this) == rhs );}
        };
        
        //! Generate the coil geom
        bool run(FieldHandle& mesh, Args& args); 
    };



    class SCISHARE ModelTMSCoilDipoleAlgo : public AlgoBase
    {
      public:
        ModelTMSCoilDipoleAlgo()
        {}

        struct Args
        {
            double totalCurrent;
            size_t numberSegments;
            double coilRadiusInner;
            double coilRadiusOuter;
            double coilDistanceOuter;
            size_t coilLevelDetails;
            int type;
            
            inline bool operator==(const Args& rhs)
            { 
                return 
                totalCurrent == rhs.totalCurrent && 
                numberSegments == rhs.numberSegments &&
                coilRadiusInner == rhs.coilRadiusInner && 
                coilRadiusOuter == rhs.coilRadiusOuter && 
                coilLevelDetails == rhs.coilLevelDetails &&
                coilDistanceOuter == rhs.coilDistanceOuter &&
                type == rhs.type;
            }
            inline bool operator!=( const Args& rhs){return !( (*this) == rhs );}
        };
        
        //! Generate the coil geom
        bool run(FieldHandle& mesh, Args& args); 
    };


} // end namespace

#endif
