//Author: p.petrov 2015 

#include <Core/Algorithms/Math/TMS/ModelGenericCoilAlgo.h>

#include <Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>

#include <Core/Datatypes/FieldInformation.h>

#include <Core/Math/MiscMath.h>

#include <string>
#include <iostream>
#include <cstdlib>

#include <boost/date_time/posix_time/posix_time.hpp>

using namespace SCIRun;

//! struct holding argumetns from the terminal
struct mainargs
{
	int type;
	int lod;
	double istep;
	double radius;
};

void CreteDomain(FieldHandle& fh, const std::vector<Vector>& points)
{
	FieldInformation fi("PointCloudMesh",1,"vector");
	fi.make_pointcloudmesh();
	fi.make_lineardata();
	fi.make_vector();

	fh = CreateField(fi);

	VMesh* mesh = fh->vmesh();

	for(size_t i = 0; i < points.size(); i++)
	{
		const Point p(points[i]);
		mesh->add_point(p);
	}
}

void CreateCoil(FieldHandle& fh, int type, int lod, double rad)
{
	if(type == 1)
	{
		//! single circle coil
		SCIRunAlgo::ModelTMSCoilSingleAlgo coilAlgoSingle;

		SCIRunAlgo::ModelTMSCoilSingleAlgo::Args coilAlgoSingleArgs;
		coilAlgoSingleArgs.type = 1; //single figure-0 shaped coil
		coilAlgoSingleArgs.wireCurrent = 55.0d;//55 A
		coilAlgoSingleArgs.coilRadius = rad;//0.044d;
		coilAlgoSingleArgs.coilDistanceOuter = 0.002d;
		coilAlgoSingleArgs.coilLevelDetails = lod;

		coilAlgoSingle.run(fh,coilAlgoSingleArgs);
	}

	if(type == 2)
	{
		//! spiral coil
		SCIRunAlgo::ModelTMSCoilSpiralAlgo coilAlgoSpiral;

		SCIRunAlgo::ModelTMSCoilSpiralAlgo::Args coilAlgoSpiralArgs;
		coilAlgoSpiralArgs.type = 1; //single figure-0 shaped coil
		coilAlgoSpiralArgs.wireCurrent = 55.0d;//55 A
		coilAlgoSpiralArgs.wireLoops = 9; //????????????????????????????????????????????????????????
		coilAlgoSpiralArgs.coilRadiusInner = 0.029d;
		coilAlgoSpiralArgs.coilRadiusOuter = 0.044d;
		coilAlgoSpiralArgs.coilDistanceOuter = 0.002d;
		coilAlgoSpiralArgs.coilLevelDetails = lod;

		coilAlgoSpiral.run(fh,coilAlgoSpiralArgs);
	}

	if(type == 3)
	{
		//! dipoles coil
		SCIRunAlgo::ModelTMSCoilDipoleAlgo coilAlgoDipole;
		
		SCIRunAlgo::ModelTMSCoilDipoleAlgo::Args coilAlgoDipoleArgs;
		coilAlgoDipoleArgs.type = 1; //single figure-0 shaped coil
		coilAlgoDipoleArgs.totalCurrent = 55.0d;//55 A
		coilAlgoDipoleArgs.numberSegments = 16; //????????????????????????????????????????????????????????
		coilAlgoDipoleArgs.coilRadiusInner = 0.029d;
		coilAlgoDipoleArgs.coilRadiusOuter = 0.044d;
		coilAlgoDipoleArgs.coilDistanceOuter = 0.002d;
		coilAlgoDipoleArgs.coilLevelDetails = lod;

		coilAlgoDipole.run(fh,coilAlgoDipoleArgs);
	}
}

void CalcAnalyticalBfield(const std::vector<Vector>& points, std::vector<Vector>& data, const double R,const double current)
{
	/// the analytical B-field: B = MU_0*I / 2 * (  R^2 / ( X^2+ R^2 ) ^ 3/2 )
	/// where R is the radius of the circular coil and X is an offset along its central/middle axis
	for(size_t i =0; i < points.size(); i++)
	{
		Vector v;
		double Z = points[i].z();
		v.z( (current * 2.0 * M_PI * 1E-7) * ( R*R / Pow( (Z*Z + R*R), 3.0 / 2.0 ) ) );
		data.push_back(v);
	}

}

void PrintAnalyticalBfield(const std::vector<Vector>& data)
{

	for(size_t i =0; i < data.size(); i++)
	{
		std::cout << data[i].length() << ", ";
	}

	std::cout << std::endl;
}

void printUsageInfo(char *progName) 
{
	std::cerr << "\n Usage: " << progName << " type lod istep radius\n\n";
	std::cerr << "\t Reports magnetic field magnitude at axial distance of 1, 2, 3, 5, 8 cm.  \n";
	std::cerr << "\t where, type: 1 - circular coil  \n";
	std::cerr << "\t              2 - spiral coil \n";
	std::cerr << "\t              3 - dipoles model \n";
	std::cerr << "\t              0 - analytical \n";
	std::cerr << "\t lod: level of detail for coil geom [integer] \n";
	std::cerr << "\t istep: integration step [real]\n";
	std::cerr << "\t radius: radius of the coil [real] \n\n";
}

bool isNumber(const std::string& s)
{
	return !s.empty() && s.find_first_not_of("-.0123456789") == std::string::npos;
}

int parseArgs(mainargs& margs, int argc, char *argv[]) 
{
	std::string arg1(argv[1]);
	std::string arg2(argv[2]);
	std::string arg3(argv[3]);
	std::string arg4(argv[4]);

	if( !isNumber(arg1) || !isNumber(arg2) || !isNumber(arg3) || !isNumber(arg4))
	{
		return 0;
	}

	//margs.type = std::stoi(arg1);
	//margs.lod = std::stoi(arg2);
	//margs.istep = std::stod(arg3);
	margs.type = atoi(arg1.c_str());
	margs.lod = atoi(arg2.c_str());
	margs.istep = atof(arg3.c_str());
	margs.radius = atof(arg4.c_str());

	//check valid range for type
	if(margs.type < 0 || margs.type > 3)
	{
		return 0;
	}

	return 1;
}

void PrintResults(MatrixHandle& mh, long timing)
{

	const DenseMatrix* const m = dynamic_cast<DenseMatrix *>(mh.get_rep());

	//size_type s = m->get_data_size();
	Vector v;

	for(index_type row = 0; row < m->nrows(); row++)
	{
		v.x(m->get(row,0));
		v.y(m->get(row,1));
		v.z(m->get(row,2));
		
		//std::cout << "Vector: " << v << " len=" << v.length() << std::endl;
		std::cout << v.length() << ", ";
	}

	std::cout << timing ;

	std::cout << std::endl;
}

int main(int argc, char *argv[])
{
	mainargs margs;

	if (argc < 5 || argc > 5) 
	{
		printUsageInfo(argv[0]);
		return 2;
	}

	if (!parseArgs(margs, argc, argv)) 
	{
		printUsageInfo(argv[0]);
		return 2;
	}

	//!sample points are fixed
	std::vector<Vector> domainPoints;
	domainPoints.reserve(5);

	Vector p1(0.0d, 0.0d, 0.01d);
	Vector p2(0.0d, 0.0d, 0.02d);
	Vector p3(0.0d, 0.0d, 0.03d);
	Vector p4(0.0d, 0.0d, 0.05d);
	Vector p5(0.0d, 0.0d, 0.08d);
	domainPoints.push_back(p1);
	domainPoints.push_back(p2);
	domainPoints.push_back(p3);
	domainPoints.push_back(p4);
	domainPoints.push_back(p5);

	if(margs.type == 0)
	{
		std::vector<Vector> dataPoints;
		dataPoints.reserve(5);

		CalcAnalyticalBfield(domainPoints, dataPoints, margs.radius, 55);

		PrintAnalyticalBfield(dataPoints);

		return 0;
	}

	//! timing the process: start
	boost::posix_time::ptime pt1 = boost::posix_time::microsec_clock::local_time();


	//! create the domain
	FieldHandle domainMesh;
	MatrixHandle domainData;
	CreteDomain(domainMesh,domainPoints);

	//! creates coil
	FieldHandle coilMesh;
	CreateCoil(coilMesh, margs.type, margs.lod, margs.radius);

	//! the generic solver for the magnetic field
	SCIRunAlgo::BiotSavartSolverAlgo solverAlgo;

	//! experiment with explicitly definig integration step
	solverAlgo.SetIntegrationStep(margs.istep);

	//! solver 
	solverAlgo.run(domainMesh, coilMesh, 1, domainData);

	//! timing the process: end
	boost::posix_time::ptime pt2 = boost::posix_time::microsec_clock::local_time();

	//! timing the duration im milliseconds
	boost::posix_time::time_duration pt_diff = pt2 - pt1;
	
	//! export results
	PrintResults( domainData, pt_diff.total_milliseconds() );
	//std::string  strdata = to_string(domainData);
	//std::cout << strdata << std::endl;

	return 0;
}
