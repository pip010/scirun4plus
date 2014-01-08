/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


/*
 *  StreamlineAnalyzer.h:
 *
 *  Written by:
 *   Allen Sanderson
 *   SCI Institute
 *   University of Utah
 *   September 2006
 *
 */

#if !defined(StreamlineAnalyzerLib_h)
#define StreamlineAnalyzerLib_h

#include <iostream>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

#include <vector>

// Basic interface between the outside world and the base libs.
void getFieldlineProperties( double pts[][3],
			     unsigned int npts,
			     unsigned int toroidalWinding,
			     unsigned int poloidalWinding,
			     unsigned int islands );


enum FieldlineType { UNKNOWN  = 0,

		     PERIODIC = 10,
		     RATIONAL = 11,
		     O_POINT  = 12,
		     X_POINT  = 13,
		     
		     QUASI_PERIODIC = 20,
		     IRRATIONAL     = 21,
		     ISLAND_CHAIN   = 22,
		     
		     CHAOTIC = 30 };
  
struct FieldlineInfo {
  FieldlineType type;

  unsigned int toroidalWinding;
  unsigned int poloidalWinding;
  unsigned int skip;
  unsigned int islands;
  float nnodes;

  float confidence;

  unsigned int zIndexOffset;
};

namespace FusionPSE {

using namespace SCIRun;

class FieldlineLib
{
public:

  Point interpert( Point lastPt, Point currPt, double t );

  int ccw( Vector v0, Vector v1 );

  int intersect( Point l0_p0, Point l0_p1,
		 Point l1_p0, Point l1_p1 );

  void convexHull( vector< Point > &hullPts,
		   vector< unsigned int > &ordering,
		   unsigned int &m,
		   unsigned int toroidalWinding,
		   int dir );

  bool hullCheck( vector< Point > &points,
		  unsigned int toroidalWinding );

  unsigned int factorial( unsigned int n0, unsigned int n1 );

  Point circle(Point &pt1, Point &pt2, Point &pt3);

  bool IsPerpendicular(Point &pt1, Point &pt2, Point &pt3);

  Point CalcCircle(Point &pt1, Point &pt2, Point &pt3);

  bool
  IntersectCheck( vector< Point >& points, unsigned int nbins );

  unsigned int Blankinship( unsigned int toroidalWinding,
			    unsigned int poloidalWinding,
			    unsigned int offset = 1 );

  void poloidalWindingCheck( vector< Point > &points,
		   vector< unsigned int > &poloidalWindingset,
		   unsigned int maxToroidalToroidalWinding,
		   unsigned int &bestToroidalWinding,
		   unsigned int &bestPoloidalWinding,
		   double &bestHitrate,
		   unsigned int &nextBestToroidalWinding,
		   unsigned int &nextBestPoloidalWinding,
		   double &nextBestHitrate,
		   unsigned int level );

  double
  poloidalWindingStats( vector< Point >& poloidalWinding_points,
	      unsigned int poloidalWinding );

  bool
  rationalCheck( vector< Point >& points,
		 unsigned int toroidalWinding,
		 unsigned int &island,
		 float &avenode,
		 float delta=0.01);

  bool
  islandChecks( vector< Point >& points,
		unsigned int toroidalWinding,
		unsigned int &islands,
		float &avenode );

  bool
  basicChecks( vector< Point >& points,
	       Vector & globalCentroid,
	       unsigned int &toroidalWinding,
	       unsigned int &poloidalWinding,
	       unsigned int &skip,
	       unsigned int &type,
	       unsigned int &islands,
	       float &avenode,
	       bool &groupCCW,
	       unsigned int &toroidalWindingNextBest );


  FieldlineInfo 
  fieldlineProperties( vector< Point > &ptList,
		       unsigned int override,
		       unsigned int maxToroidalWinding,
		       float hitrate );

  FieldlineInfo 
  fieldlineProperties( vector< Point > &points,
		       unsigned int maxToroidalToroidalWinding );


  unsigned int
  islandProperties( vector< Point > &points,
		    Vector &baseCentroid,
		    unsigned int &startIndex,
		    unsigned int &middleIndex,
		    unsigned int &stopIndex,
		    unsigned int &nodes );

  unsigned int
  surfaceOverlapCheck( vector< vector< Point > > &bins,
		    unsigned int toroidalWinding,
		    unsigned int skip,
		    unsigned int &nnodes );

  unsigned int
  surfaceGroupCheck( vector< vector< Point > > &bins,
		     unsigned int i,
		     unsigned int j,
		     unsigned int nnodes );

  unsigned int
  removeOverlap( vector< vector < Point > > &bins,
		 unsigned int &nnodes,
		 unsigned int toroidalWinding,
		 unsigned int poloidalWinding,
		 unsigned int skip,
		 unsigned int island );

  unsigned int
  smoothCurve( vector< vector < Point > > &bins,
	       unsigned int &nnodes,
	       unsigned int toroidalWinding,
	       unsigned int poloidalWinding,
	       unsigned int skip,
	       unsigned int island );

  unsigned int
  mergeOverlap( vector< vector < Point > > &bins,
		unsigned int &nnodes,
		unsigned int toroidalWinding,
		unsigned int poloidalWinding,
		unsigned int skip,
		unsigned int island );

template<class T>
void 
standardDeviationCheck( vector< T > &sampleSet,
			unsigned int maxOffset,
			unsigned int &besOffset,
			double &bestSD,
			unsigned int sampleType );
};


template<class T>
void FieldlineLib::
standardDeviationCheck( vector< T > &sampleSet,
			unsigned int maxOffset,
			unsigned int &besOffset,
			double &bestSD,
			unsigned int sampleType )
{
  besOffset = 0;
  bestSD = 1.0e9;

  unsigned int nSamples = sampleSet.size();

  // The premise is that for a given poloidal winding the number of
  // points should be consistant between each Nth punction point,
  // where N is the poloidal winding. For instance, if the poloidal
  // winding is 5 and the number is 20. Then the pattern
  // should be:

  // 0 1 1 1 2 - 2 3 3 3 4 - 4 5 5 5 6

  // In this case the different between every 5th value (the poloidal
  // winding) should be 2 (the poloidal winding).

  for( unsigned int offset=1;
       (offset<maxOffset && offset<nSamples/2);
       ++offset )
  {
    double average_sd = 0;
    
    for( unsigned int j=0; j<offset; ++j )
    {
      // Get the average value of samples between each offset.
      double sampleAve = 0;
      unsigned int npts = 0;
      
      double sumofsquares = 0;

      if( sampleType == 0 )
      {
	for( unsigned int i=j; i<nSamples-offset; i+=offset)
        {
	  sampleAve += (sampleSet[i+offset] - sampleSet[i]);
	  ++npts;
	}
	
	sampleAve /= (float) npts;
	
	// Get the standard deviation for the samples
	for( unsigned int i=j; i<nSamples-offset; i+=offset)
	{
	  double diff = (sampleSet[i+offset] - sampleSet[i]) - sampleAve;
	  
	  sumofsquares += diff * diff;
	}
      }
      else // if( sampleType == 1 )
      {
	for( unsigned int i=j; i<nSamples; i+=offset)
        {
	  sampleAve += sampleSet[i];
	  ++npts;
	}
	
	sampleAve = sampleAve / (float) npts;
	
	// Get the standard deviation for the samples
	for( unsigned int i=j; i<nSamples; i+=offset)
	{
	  double diff = sampleSet[i] - sampleAve;
	  
	  sumofsquares += (diff * diff);
	}
      }

      // Get the standard deviation for this grouping.
      average_sd += sqrt( sumofsquares / (float) npts );
    }

    average_sd = average_sd / (float) offset;

    if( bestSD > average_sd ) {

      bestSD = average_sd;

      besOffset = offset;
    }
  }
}


} // End namespace FusionPSE

#endif // StreamlineAnalyzerLib_h
