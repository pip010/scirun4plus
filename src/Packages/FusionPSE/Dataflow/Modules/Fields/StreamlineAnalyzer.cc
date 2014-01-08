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
 *  StreamlineAnalyzer.cc:
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   September 2005
 *
 */

#include <Core/Thread/Mutex.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Packages/FusionPSE/Dataflow/Modules/Fields/StreamlineAnalyzerLib.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Core/Containers/Handle.h>

namespace FusionPSE {

using namespace SCIRun;

class StreamlineAnalyzer : public Module
{
public:
  StreamlineAnalyzer(GuiContext* ctx);
  virtual ~StreamlineAnalyzer();

  virtual void execute();

protected:
  GuiString gui_PlanesStr_;
  GuiInt gui_PlanesInt_;
  GuiInt gui_Color_;
  GuiInt gui_MaxToroidalWinding_;
  GuiInt gui_Override_;
  GuiDouble gui_HitRate_;
  GuiInt gui_CurveMesh_;
  GuiInt gui_AdjustPlane_;
  GuiInt gui_ShowIslands_;
  GuiInt gui_Overlaps_;

  vector< double > planes_;
};

class StreamlineAnalyzerAlgo
{
  public:
    Mutex lock;
    
    StreamlineAnalyzerAlgo() :
      lock("StreamlineAnalyzerAlgo") {}
  
    //! virtual interface. 
    void execute(FieldHandle& src,
		 FieldHandle& dst,
		 FieldHandle& pccdst,
		 FieldHandle& pcsdst,
		 FieldHandle& tmpdst,
		 vector< double > &planes,
		 unsigned int color,
		 unsigned int showIslands,
		 unsigned int overlaps,
		 unsigned int maxToroidalWinding,
		 unsigned int override,
		 float hitrate,
		 vector< pair< unsigned int,
		 unsigned int > > &topology,
		 bool is_curvemesh,
		 int adjust_plane );

  virtual ~StreamlineAnalyzerAlgo() {}

  virtual void
  loadCurve( FieldHandle &field_h,
	     vector< vector < vector < Point > > > &nodes,
	     unsigned int color,
	     double color_value );

  virtual void
  loadSurface( FieldHandle &field_h,
	       vector < vector < vector < Point > > > &nodes,
	       unsigned int nnodes,
	       unsigned int islands,
	       unsigned int skip,
	       unsigned int adjust_plane,
	       unsigned int color,
	       double color_value );

  short SIGN( double val ) { return (val < 0 ? -1 : 1 ); };
           
  protected:
    FieldlineLib FLlib;         

    unsigned int safetyFactorConsistant;
    unsigned int poloidalPeriodicyMatch;
    unsigned int poloidalPeriodicyMismatch;
};


DECLARE_MAKER(StreamlineAnalyzer)

StreamlineAnalyzer::StreamlineAnalyzer(GuiContext* context)
  : Module("StreamlineAnalyzer", context, Source, "Fields", "FusionPSE"),
    gui_PlanesStr_(context->subVar("planes-list"), "0.0"),
    gui_PlanesInt_(context->subVar("planes-quantity"), 0),
    gui_Color_(context->subVar("color"), 1),
    gui_MaxToroidalWinding_(context->subVar("maxToroidalWinding"), 30),
    gui_Override_(context->subVar("override"), 0),
    gui_HitRate_(context->subVar("hitrate"), 0.90),
    gui_CurveMesh_(context->subVar("curve-mesh"), 1),
    gui_AdjustPlane_(context->subVar("plane-adjust"), -1),
//    gui_ScalarField_(context->subVar("scalar-field"), 1),
    gui_ShowIslands_(context->subVar("show-islands"), 0),
    gui_Overlaps_(context->subVar("overlaps"), 1)
{
}

StreamlineAnalyzer::~StreamlineAnalyzer()
{
}

void
StreamlineAnalyzer::execute()
{
  // The streamline field input is required.
  FieldHandle sl_field_input_handle;

  get_input_handle( "Input Streamlines", sl_field_input_handle, true );
  FieldInformation sl_fi(sl_field_input_handle);
  
  if (!(sl_fi.is_curvemesh()))
  { 
    error("Only available for CurveFields.");
    return;
  }

  if (!(sl_fi.is_scalar()))
  {
    error("Only available for Scalar data.");
    return;
  }

  // Get the numbers of planes to display the resulting Poincare plot.
  std::vector< double > planes(0);

  if( gui_PlanesInt_.get()) 
  {
    unsigned int nplanes = gui_PlanesInt_.get();

    for( unsigned int i=0; i<nplanes; i++ )
      planes.push_back(2.0 * M_PI * (double) i / (double) nplanes );
  } 
  else 
  {
    std::istringstream plist(gui_PlanesStr_.get());
    double plane;
    while(!plist.eof()) 
    {
      plist >> plane;
      if (plist.fail()) 
      {
        if (!plist.eof()) 
        {
          plist.clear();
          warning("List of Planes was bad at character " +
            to_string((int)(plist.tellg())) +
            "('" + ((char)(plist.peek())) + "').");
        }
        break;
      } 
      else if (!plist.eof() && plist.peek() == '%') 
      {
        plist.get();
        plane = 0 + (2.0*M_PI - 0) * plane / 100.0;
      }

      if( 0 <= plane && plane <= 2.0*M_PI )
        planes.push_back(plane);
      else 
      {
        error("Plane is not in the range of 0 to 2 PI.");
        return;
      }
    }
  }

  // Check the plane list for differences.
  if( planes_.size() != planes.size() )
  {
    inputs_changed_ = true;

    planes_.resize(planes.size());

    for( unsigned int i=0; i<planes.size(); ++i )
      planes_[i] = planes[i];

  } 
  else 
  {
    for( unsigned int i=0; i<planes.size(); ++i )
    {
      if( fabs( planes_[i] - planes[i] ) > 1.0e-4 ) 
      {
        planes_[i] = planes[i];
        inputs_changed_ = true;
      }
    }
  }

  std::cerr << "StreamlineAnalyzer executing " << std::endl;

  // If no data or a changed recalcute.
  if( inputs_changed_ ||

      !oport_cached("Output Poincare") ||

      gui_Color_.changed( true ) ||
      gui_MaxToroidalWinding_.changed( true ) ||
      gui_Override_.changed( true ) ||
      gui_HitRate_.changed( true ) ||
      gui_CurveMesh_.changed( true ) ||
      gui_AdjustPlane_.changed( true ) ||
      gui_ShowIslands_.changed( true ) ||
      gui_Overlaps_.changed( true ) ) 
  {
    update_state( Executing );

    std::vector< std::pair< unsigned int, unsigned int > > topology;

    FieldHandle sl_field_output_handle;
    FieldHandle pcc_field_output_handle;
    FieldHandle pcs_field_output_handle;
    FieldHandle pts_field_output_handle;

    StreamlineAnalyzerAlgo algo;

    algo.execute(sl_field_input_handle,   sl_field_output_handle,
		 pcc_field_output_handle, pcs_field_output_handle,
		 pts_field_output_handle,
		 planes_,
		 gui_Color_.get(),
		 gui_ShowIslands_.get(),
		 gui_Overlaps_.get(),
		 gui_MaxToroidalWinding_.get(), gui_Override_.get(),
		 gui_HitRate_.get(),
		 topology,
		 gui_CurveMesh_.get(),
		 gui_AdjustPlane_.get());

    // Send the data downstream
    send_output_handle( "Output Poincare",     sl_field_output_handle );
    send_output_handle( "Output Centroids",    pcc_field_output_handle );
    send_output_handle( "Output Separatrices", pcs_field_output_handle );
    send_output_handle( "Output Points",       pts_field_output_handle );

    std::cerr << "StreamlineAnalyzer done " << std::endl;
  }
}


void
StreamlineAnalyzerAlgo::
execute(FieldHandle& ifield_h,
	FieldHandle& ofield_h,
	FieldHandle& opccfield_h,
	FieldHandle& opcsfield_h,
	FieldHandle& optsfield_h,
	std::vector< double > &planes,
	unsigned int color,
	unsigned int showIslands,
	unsigned int overlaps,
	unsigned int maxToroidalWinding,
	unsigned int override,
	float hitrate,
	std::vector< std::pair< unsigned int, unsigned int > > &topology,
	bool is_curvemesh,
	int adjust_plane)
{
  VField *ifield = ifield_h->vfield();
  VMesh  *imesh =  ifield->vmesh();

  FieldInformation fo(ifield_h);
  fo.make_double();
  
  // Make the field now as it unstructured (and irregular).
  if (is_curvemesh)
  {
    fo.make_curvemesh();
  }
  else
  {
    fo.make_quadsurfmesh();
  }

  ofield_h = CreateField(fo);

  // Create two point cloud meshes for the centroids and separatices.
  FieldInformation fi_pc("PointCloudMesh",0,"double");

  // Create the centroids.
  opccfield_h = CreateField(fi_pc);
  
  VMesh  *opccmesh = opccfield_h->vmesh();
  VField *opccfield = opccfield_h->vfield();

  // Create the separatices.
  opcsfield_h = CreateField(fi_pc);
  
  VMesh  *opcsmesh = opcsfield_h->vmesh();
  VField *opcsfield = opcsfield_h->vfield();

  // Create the points.
  fo.make_curvemesh();
  optsfield_h = CreateField(fo);
  
  VMesh  *optsmesh = optsfield_h->vmesh();
  VField *optsfield = optsfield_h->vfield();

  VMesh::Node::array_type enodes(2);
  VMesh::Node::index_type n1, n2;

  // Input iterators
  VMesh::Node::iterator inodeItr, inodeEnd;
  VMesh::Node::index_type inodeNext;
  std::vector< VMesh::Node::index_type > inodeGlobalStart;

  imesh->begin( inodeItr );
  imesh->end( inodeEnd );

  if(inodeItr == inodeEnd) return;

  std::vector< std::vector< Point > > fieldlines;
  std::vector< Point > *ptList = new std::vector< Point >;

  // The index of the first node for the next streamline.
  std::ostringstream str;
  str << "Streamline " << fieldlines.size()+1 << " Node Index";

  unsigned int index;
      
  if( ifield->get_property( str.str(), index ) )
    inodeNext = index;
  else
    inodeNext = *inodeEnd;

  // Spearate each fieldline into a simple point list.
  while (inodeItr != inodeEnd) 
  {
    // Get the next point.
    Point pt;
    imesh->get_center(pt, *inodeItr);
    ptList->push_back(pt);

    ++inodeItr;

    // Next point is the start of the next fieldline.
    if( *inodeItr == inodeNext)
    {
      // Store this fieldline point list.
      if( ptList->size() )
      {
	fieldlines.push_back( *ptList );
      }

      if( inodeItr != inodeEnd ) 
      {
	// The index of the first node for the next streamline.
	std::ostringstream str;
	str << "Streamline " << fieldlines.size()+1 << " Node Index";
	
	if( ifield->get_property( str.str(), index ) )
	  inodeNext = index;
	else
	  inodeNext = *inodeEnd;
	
	ptList = new std::vector< Point >;
      }
    }
  }

  safetyFactorConsistant = 0;
  poloidalPeriodicyMatch = 0;
  poloidalPeriodicyMismatch = 0;


  // Now analyze and bin the points for each fieldline.
  for( unsigned int i=0; i<fieldlines.size(); i++ ) 
  {
    std::cerr << " STARTING STREAMLINE " << i << std::endl << std::endl;
    
    // Get the winding information for each fieldline.
    FieldlineInfo fieldlineInfo =
      FLlib.fieldlineProperties( fieldlines[i],
				 override,
				 maxToroidalWinding,
				 hitrate );
    
    FieldlineType type           = fieldlineInfo.type;
    unsigned int toroidalWinding = fieldlineInfo.toroidalWinding;
    unsigned int poloidalWinding = fieldlineInfo.poloidalWinding;
    unsigned int islands         = fieldlineInfo.islands;
    unsigned int skip            = fieldlineInfo.skip;
    unsigned int nnodes          = (unsigned int) fieldlineInfo.nnodes;
    unsigned int zIndexOffset    = fieldlineInfo.zIndexOffset;

    bool completeIslands = true;

    if( toroidalWinding == 0 ) 
    {
      std::cerr << i << " SKIPPING TOROIDALWINDING OF 0" << std::endl;

      std::pair< unsigned int, unsigned int > topo( 0, 0 );
      topology.push_back(topo);

      continue;
    }

    else if( type == UNKNOWN || type == CHAOTIC ) 
    {
      std::cerr << i << " SKIPPING " << std::endl;

      std::pair< unsigned int, unsigned int > topo( 0, 0 );
      topology.push_back(topo);

      continue;
    }

    // Get the direction of the streamline toroidalWinding.
    Point lastPt = fieldlines[i][0];
    Point currPt = fieldlines[i][1];

    bool CCWstreamline = (atan2( lastPt.y(), lastPt.x() ) <
			  atan2( currPt.y(), currPt.x() ));


    // Put all of the points into the bins for each plane.
    std::vector< std::vector< std::vector < Point > > > puncturePts;

    puncturePts.resize( planes.size() );

    unsigned int startIndex = 0;

    for( unsigned int p=0; p<planes.size(); ++p ) 
    {
      Vector planeN;
      Vector planePt(0,0,0);

      // Go through the planes in the same direction as the streamline.
      if( CCWstreamline )
      {
	planeN = Vector( cos(planes[p]),
			 sin(planes[p]),
			 0 );
      }
      else
      {
	planeN = Vector( cos(planes[planes.size()-1-p]),
			 sin(planes[planes.size()-1-p]),
			 0 );
      }

      // Set up the plane equation.
      double plane[4];
      
      plane[0] = planeN.x();
      plane[1] = planeN.y();
      plane[2] = planeN.z();
      plane[3] = Dot( planePt, planeN);

//      cerr << "Plane " << p << " is " << plane << endl;
        
//    cerr << "Starting new streamline binning " << c << endl;

      puncturePts[p].resize( toroidalWinding );
      int bin = 0;

      // So to get the winding groups consistant start examining the
      // streamline in the same place for each plane.
      Vector lastPt, currPt( fieldlines[i][startIndex] );
      double lastDist, currDist = Dot( planeN, currPt ) - plane[3];


      // Set up the Z plane equation.
      Vector planeNZ( 0, 0, 1 );

      double planeZ[4];
      
      planeZ[0] = planeNZ.x();
      planeZ[1] = planeNZ.y();
      planeZ[2] = planeNZ.z();
      planeZ[3] = Dot(planePt, planeNZ);

      double lastDistZ, currDistZ = Dot(planeNZ, currPt) - planeZ[3];

      unsigned jj = 0;

      for( unsigned int j=startIndex+1; j<fieldlines[i].size(); ++j )
      {
	lastPt = currPt;
	currPt = Vector(fieldlines[i][j]);

	lastDist = currDist;
	currDist = Dot( planeN, currPt ) - plane[3];

	// First look at only points that intersect the plane.
	if( SIGN(lastDist) != SIGN(currDist) ) 
	{
	  Vector dir(currPt-lastPt);

	  double dot = Dot(planeN, dir);

	  // If the segment is in the same direction as the plane then
	  // find where it intersects the plane.
	  if( dot > 0.0 )
	  {
	    // So to get the winding groups consistant start examining
	    // the streamline in the same place for each plane so
	    // store the index of the first puncture point.
	    if( p == 0 && puncturePts[p][bin].size() == 0 )
	      startIndex = j - 1;

	    Vector w = lastPt - planePt;

	    double t = -Dot(planeN, w ) / dot;

	    Point point = Point(lastPt + dir * t);
	  
	    puncturePts[p][bin].push_back( point );
	    
	    bin = (bin + 1) % toroidalWinding;
	  }
	}

	// Find the positive zero crossings which indicate a poloidal winding.
	if( p == 0 && puncturePts[p][bin].size() &&
	    j+zIndexOffset < fieldlines[i].size() )
	{
	  // Poloidal plane distances.
	  lastDistZ = currDistZ;
	  currDistZ = Dot( planeNZ, currPt ) - planeZ[3];
	  
	  // First look at only points that intersect the toroiadal plane.
	  if( SIGN(lastDistZ) != SIGN(currDistZ) ) 
	  {
	    Vector dir(currPt-lastPt);
	    
	    double dot = Dot(planeNZ, dir);
	    
	    // If the segment is in the same direction as the toroidal plane
	    // then find where it intersects the plane.
	    if( dot > 0.0 )
	    {
	      Vector w = (Vector) lastPt - planePt;
	      
	      double t = -Dot(planeNZ, w ) / dot;
        
	      Point point = Point(lastPt + dir * t);
	      
	      n1 = n2;
	      n2 = optsmesh->add_node(fieldlines[i][j+zIndexOffset]);

	      optsfield->resize_fdata();
	      optsfield->set_value( n2, n2 );
	      
	      jj = j;
	      
	      // Add the next element (i.e. an edge)
	      if( optsmesh->num_nodes() > 0 )
	      {
		enodes[0] = n1;
		enodes[1] = n2;
		optsmesh->add_elem(enodes);
	      }
	    }
	  }
	}
      }
    }
    
    bool VALID = true;

    // Sanity check
    for( unsigned int p=0; p<planes.size(); ++p ) 
    {
      for( unsigned int j=0; j<toroidalWinding; ++j ) 
      {
        if( nnodes > puncturePts[p][j].size() )
          nnodes = puncturePts[p][j].size();

        if( puncturePts[p][j].size() < 1 ) 
        {
          std::cerr << "INVALID - Plane " << p
               << " bin  " << j
               << " number of points " << puncturePts[p][j].size()
               << std::endl;
          VALID = false;

          return;
        }
      }
    }

    // Get the rest of the info only from the phi = zero plane.
    unsigned int p;

    if( CCWstreamline )
      p = 0;
    else
      p = planes.size()-1;

    // Get the centroid of each toroidal winding group and all
    // puncture points.
//    Vector globalCentroid(0,0,0);
    std::vector< Vector > localCentroids;
    std::vector< Vector > localSeparatrices[2];

    localCentroids.resize(toroidalWinding);
    localSeparatrices[0].resize(toroidalWinding);
    localSeparatrices[1].resize(toroidalWinding);

    for( unsigned int j=0; j<toroidalWinding; ++j ) 
    {
      localCentroids[j] = Vector(0,0,0);
      
      for( unsigned int k=0; k<puncturePts[p][j].size(); ++k ) 
        localCentroids[j] += (Vector) puncturePts[p][j][k];

      if( puncturePts[p][j].size() ) 
      {
        localCentroids[j] /= (double) puncturePts[p][j].size();

//        globalCentroid += localCentroids[j];
      }
    }

//    globalCentroid /= toroidalWinding;

    // Get the direction of the points within a group.
//    Vector v0 = (Vector) puncturePts[p][0][0] - globalCentroid;
//    Vector v1 = (Vector) puncturePts[p][0][1] - globalCentroid;

//    bool groupCCW = (FLlib.ccw( v0, v1 ) == 1);
//    cerr << 0.0<< "  " << groupCCW << endl;

    if( type == ISLAND_CHAIN ) 
    {
      for( unsigned int j=0; j<toroidalWinding; ++j ) 
      {
        unsigned int startIndex;
        unsigned int middleIndex;
        unsigned int stopIndex;
        unsigned int nodes;
        
        Vector localCentroid;

        unsigned int turns =
          FLlib.islandProperties( puncturePts[p][j], localCentroid,
              startIndex, middleIndex, stopIndex, nodes );
	
// 	cerr << "Island " << i  << "   "
// 	     << "Turns " << turns  << "   "
// 	     << "nodes " << nodes  << "   "
// 	     << "Indexes "
// 	     << startIndex  << "  "
// 	     << middleIndex << "  "
// 	     << stopIndex   << endl;


//	if( turns < 3 )
//	completeIslands = false;

        if( turns >= 2 ) 
        {
// 	  localSeparatrices[0][j] = (Vector) puncturePts[p][j][startIndex];
// 	  localSeparatrices[1][j] = (Vector) puncturePts[p][j][middleIndex];
        }

        if( turns == 3 ) 
        {
          unsigned int index0 = (middleIndex - startIndex ) / 2;
          unsigned int index1 = (  stopIndex - middleIndex) / 2;

          // 	  cerr << "Indexes mid " << nodes << " nodes "
          // 	       << "  " << ( startIndex + index0)%nodes 
          // 	       << "  " << (middleIndex - index0)%nodes
          // 	       << "  " << (middleIndex + index1)%nodes
          // 	       << "  " << (  stopIndex - index1)%nodes << endl;

          localCentroids[j] =
            ( (Vector) puncturePts[p][j][( startIndex + index0)%nodes] + 
              (Vector) puncturePts[p][j][(middleIndex - index0)%nodes] +

              (Vector) puncturePts[p][j][(middleIndex + index1)%nodes] + 
              (Vector) puncturePts[p][j][(  stopIndex - index1)%nodes] ) / 4.0;
        } 
        else if( turns == 2 ) 
        {
          unsigned int index0 = (middleIndex - startIndex ) / 2;

          //  	  cerr << "Indexes mid " << nodes << " nodes "
          // 	       << "  " << ( startIndex + index0)%nodes 
          // 	       << "  " << (middleIndex - index0)%nodes
          // 	       << "  " << (middleIndex + index1)%nodes
          // 	       << "  " << (  stopIndex - index1)%nodes << endl;

          localCentroids[j] =
            ( (Vector) puncturePts[p][j][( startIndex + index0)%nodes] + 
              (Vector) puncturePts[p][j][(middleIndex - index0)%nodes] ) / 2.0;

        } 
        else if( turns == 1 ) 
        {
          unsigned int index0 = (stopIndex - startIndex ) / 2;

          //  	  cerr << "Indexes mid " << nodes << " nodes "
          // 	       << "  " << ( startIndex + index0)%nodes 
          // 	       << "  " << (middleIndex - index0)%nodes
          // 	       << "  " << (middleIndex + index1)%nodes
          // 	       << "  " << (  stopIndex - index1)%nodes << endl;

          localCentroids[j] =
            ( (Vector) puncturePts[p][j][(startIndex + index0)%nodes] + 
              (Vector) puncturePts[p][j][( stopIndex - index0)%nodes] ) / 2.0;
        }
	
// 	// Get the principal axes of the island.
// 	Vector localCentroid(0,0,0);

// 	for( unsigned int k=0; k<puncturePts[p][j].size(); ++k )
// 	  localCentroid += (Vector) puncturePts[p][j][k];

// 	localCentroid /= (float) puncturePts[p][j].size();

// 	float Ixx = 0.0;
// 	float Ixz = 0.0;
// 	float Izz = 0.0;

// 	double maxDist = 0;

// 	for( unsigned int k=0; k<puncturePts[p][j].size(); k++ ) {

// 	  Vector vec = (Vector) puncturePts[p][j][k] - localCentroid;

// 	  if( maxDist < vec.length() )
// 	    maxDist = vec.length();

// 	  Ixx += vec.z()*vec.z();
// 	  Ixz -= vec.x()*vec.z();
// 	  Izz += vec.x()*vec.x();
// 	}

// 	float alpha = atan( 2.0 * Ixz / (Ixx - Izz) ) / 2.0;

// //       cerr << "PRINCIPAL AXES " << alpha * 180.0 / M_PI << "    "
// // 	   << Ixx + Ixz * sin(alpha       )/cos(alpha       ) << "    "
// // 	   << Izz + Ixz * cos(alpha+M_PI/2)/sin(alpha+M_PI/2) << endl;

// 	if( Ixx + Ixz * sin(alpha       )/cos(alpha       ) >
// 	    Izz + Ixz * cos(alpha+M_PI/2)/sin(alpha+M_PI/2) )
// 	  localCentroid += Vector(  cos(alpha), 0, sin(alpha) ) * maxDist;
// 	else
// 	  localCentroid += Vector( -sin(alpha), 0, cos(alpha) ) * maxDist;

//	localCentroids[j] = localCentroid;

      }
    }  // if( type == ISLAND_CHAIN )

    for( unsigned int p=0; p<planes.size(); p++ ) 
    {
      if( type != RATIONAL ) 
      {
        if( overlaps == 1 || overlaps == 3 )
          FLlib.removeOverlap( puncturePts[p], nnodes,
			       toroidalWinding, poloidalWinding,
			       skip, islands );
        if( overlaps == 2 )
          FLlib.mergeOverlap( puncturePts[p], nnodes,
			      toroidalWinding, poloidalWinding,
			      skip, islands );
        else if( overlaps == 3 )
          FLlib.smoothCurve( puncturePts[p], nnodes,
			     toroidalWinding, poloidalWinding,
			     skip, islands );
      }

      bool VALID = true;

      // Sanity check
      for( unsigned int j=0; j<toroidalWinding; ++j ) 
      {
        if( nnodes > puncturePts[p][j].size() )
          nnodes = puncturePts[p][j].size();

        if( puncturePts[p][j].size() < 1 ) 
        {
          std::cerr << "INVALID - Plane " << p
               << " bin  " << j
               << " number of points " << puncturePts[p][j].size()
               << std::endl;

          VALID = false;

          return;
        }
	  
// 	cerr << "Surface " << i
// 	     << " plane " << p
// 	     << " bin " << j
// 	     << " base number of nodes " << nnodes
// 	     << " number of points " << puncturePts[p][j].size()
// 	     << endl;
      }
    }

    std::cerr << "Surface " << i << " is a " << type << "  "
       << toroidalWinding << ":" << poloidalWinding << " surface ("
       << (double) toroidalWinding / (double) poloidalWinding << ") ";
    
    if( type == ISLAND_CHAIN ) 
      std::cerr << "that contains " << islands << " islands"
        << (completeIslands ? " (Complete)" : "");
    
    std::cerr << " and has " << nnodes << " nodes" << std::endl;
        
    if( type == ISLAND_CHAIN && islands != toroidalWinding ) 
    {
      std::cerr << "WARNING - The island count does not match the toroidalWinding count" 
        << std::endl;
    }
    
    // Record the topology.
    std::pair< unsigned int, unsigned int >
      topo( toroidalWinding, poloidalWinding );

    topology.push_back(topo);


    if( islands )
    {
      if( completeIslands ) 
      {
	for( unsigned int j=0; j<toroidalWinding; ++j )
        {
	  // Centroids
	  // 	    cerr << j << "  " << localCentroids[j] << endl;
	  
	  VMesh::Node::index_type n =
	    opccmesh->add_node((Point) localCentroids[j]);
	  opccfield->resize_values();
	  opccfield->set_value( j, n );
	  
	  // Separatrices
	  unsigned int k = (j+skip) % toroidalWinding;
	  
	  unsigned int jj, kk;
	  
	  if( (localSeparatrices[0][j] - localSeparatrices[0][k]).length() <
	      (localSeparatrices[0][j] - localSeparatrices[1][k]).length() )
	    kk = 0;
	  else
	    kk = 1;

	  if( (localSeparatrices[0][j] - localSeparatrices[kk][k]).length() <
	      (localSeparatrices[1][j] - localSeparatrices[kk][k]).length() )
	    jj = 0;
	  else
	    jj = 1;

	  n = opcsmesh->add_node((Point) ((localSeparatrices[jj][j] +
					   localSeparatrices[kk][k])/2.0));

	  opcsfield->resize_values();
	  opcsfield->set_value( j, n );

//  	    n = opcsmesh->add_node((Point) localSeparatrices[0][j]);
// 	    opcsfield->resize_fdata();
// 	    opcsfield->set_value(0, n);

//  	    n = opcsmesh->add_node((Point) localSeparatrices[1][k]);
// 	    opcsfield->resize_fdata();
// 	    opcsfield->set_value(1, n);

// 	    cerr << i << "  Separatrices " << j << "  " << k << endl;
        }
      }
    } 
    else 
    { // Surface
    }

    if( !showIslands || (showIslands && islands) ) 
    {
      double color_value = 0;

      if( color == 1 )
	color_value = i;
      else if( color == 6 )
	color_value = toroidalWinding;
      else if( color == 7 )
	color_value = poloidalWinding;
      else if( color == 8 )
	color_value = (double) toroidalWinding / (double) poloidalWinding;

      // Add the points into the return field.
      if( is_curvemesh )
      {
	loadCurve( ofield_h, puncturePts, color, color_value );
      }
      else
      {
	loadSurface( ofield_h, puncturePts, nnodes, islands, skip,
		     adjust_plane, color, color_value );
      }
    }
  }

  std::cerr << std::endl << std::endl
       << "count " << fieldlines.size()
       << "  safetyFactorConsistant  " << safetyFactorConsistant
       << "  poloidalPeriodicyMatch  " << poloidalPeriodicyMatch 
       << "  poloidalPeriodicyMismatch  " << poloidalPeriodicyMismatch
       << std::endl << std::endl;
}



void
StreamlineAnalyzerAlgo::loadCurve( FieldHandle &field_h,
				   vector< vector < vector < Point > > > &nodes,
				   unsigned int color,
				   double color_value ) 
{
  lock.lock();

  VField* ofield = field_h->vfield();
  VMesh* omesh = ofield->vmesh();
  
   VMesh::Node::array_type enodes(2);
   VMesh::Node::index_type n1, n2;
    
   unsigned int nplanes = nodes.size();
   unsigned int toroidalWindings = nodes[0].size();
   
   // Loop through each plane
   for( unsigned int p=0; p<nplanes; ++p ) 
   {
     if( color == 3 )
       color_value = p;
        
     // Loop through each toroidial group
     for( unsigned int j=0; j<toroidalWindings; ++j ) 
     {
       if( color == 4 )
	 color_value = j;
       
      // Loop through each point in toroidial group
       for( unsigned int i=0; i<nodes[p][j].size(); ++i ) 
       {
	 n1 = n2;
	 n2 = omesh->add_node(nodes[p][j][i]);
	 ofield->resize_fdata();
	 
	 if( color == 2 )
	   color_value = (i*toroidalWindings+j)*nplanes + p;
	 else if( color == 5 )
	   color_value = i;
	
	 ofield->set_value( color_value, n2);

	 // Add the next element (i.e. an edge)
	 if( i > 0 )
	 {
	   enodes[0] = n1;
	   enodes[1] = n2;
	   omesh->add_elem(enodes);
	 }
       }
     }
   }   

   lock.unlock();
}
 

void
StreamlineAnalyzerAlgo::loadSurface( FieldHandle &field_h,
				     vector< vector < vector < Point > > > &nodes,
				     unsigned int nnodes,
				     unsigned int islands,
				     unsigned int skip,
				     unsigned int adjust_plane,
				     unsigned int color,
				     double color_value) 
{
  VField* ofield = field_h->vfield();
  VMesh* omesh = ofield->vmesh();

  VMesh::Node::array_type enodes(4);
  VMesh::Node::index_type n1;
  
  int base_index = -1;

  unsigned int nplanes = nodes.size();
  unsigned int toroidalWindings = nodes[0].size();
  
  int dims[2];
  
  // Add one to the first dimension to create a closed cylinder. Add
  // one to the second dimension to form a torus.
  dims[0] = nnodes + 1;
  dims[1] = nplanes * toroidalWindings + 1 - 2;
  
  // Determine if the winding group order matches the point
  // ordering. This is only needed when building surfaces.
  Vector v0 = nodes[0][   0][1] - nodes[0][0][0];
  Vector v1 = nodes[0][skip][0] - nodes[0][0][0];
  
  bool flip;
  
  if( Dot( v0, v1 ) < 0 )
    flip = true;
  else
    flip = false;
  
  // Loop through each toroidial group
  for( unsigned int j=0; j<toroidalWindings; ++j )
  {
    if( color == 4 )
      color_value = j;

    // Loop through each plane.
    for( unsigned int p=0; p<nplanes; ++p ) 
    {
      // Normally each toroidial winding group can be displayed in the
      // order received. Except for the last plane where it needs to
      // be adjusted by one group. That is if the streamline started
      // in the "correct" place. This is not always the case so it may
      // be necessary to adjust the toroidal winding group location by
      // one.
      unsigned int k, offset;

      if( p == adjust_plane )
      {
	k = (j-1 + toroidalWindings) % toroidalWindings;

	if( islands )
	  offset = 1;
	else
	  offset = 0;
      }
      else
      {
	k = j;
	offset = 0;
      }
            
      unsigned int jj = nplanes * j + p;
            
      if( color == 3 )
	color_value = jj;
            
      // Loop through each point in toroidial group.
      for(unsigned int i=0; i<nnodes; ++i )
      {
	n1 = omesh->add_node(nodes[p][k][(i+offset)%nnodes]);
	ofield->resize_fdata();

	// Each surface has a base index which is the index of the
	// first point. This is needed when multiple surfaces are
	// generted.
	if( base_index == -1 )
	  base_index = n1;

	if( color == 2 )
	  color_value = (i*toroidalWindings+j)*nplanes + p;
	else if( color == 5 )
	  color_value =  i;
	
 	ofield->set_value( color_value, n1 );

	// Add the next element (i.e. a quad).
	enodes[0] = base_index + jj     * dims[0] + i;
	enodes[1] = base_index + (jj+1) * dims[0] + i;
	enodes[2] = base_index + (jj+1) * dims[0] + i + 1;
	enodes[3] = base_index + jj     * dims[0] + i + 1;

	omesh->add_elem(enodes);
      }

      // For a surface add in the first point from the adjacent
      // toroidial group. Otherwise for an island add in the
      // first point from the current toroidal group.
      if( !islands )
      {
	if( flip )
	  k = (k-skip+toroidalWindings) % toroidalWindings;
	else
	  k = (k+skip) % toroidalWindings;
      }

      unsigned int i = nnodes;

      n1 = omesh->add_node(nodes[p][k][0]);
      ofield->resize_fdata();

      if( color == 2 )
	color_value = (i*toroidalWindings+j)*nplanes + p;
      else if( color == 5 )
	color_value =  i;
            
      ofield->set_value( color_value, n1 );
    }
  }
    
  // Add in the first toroidal group from the first plane to complete
  // the torus.
  unsigned int j = 0;
    
  if( color == 4 )
    color_value = j;
    
  // Add in the first toroidal group from the first plane to complete
  // the torus.
  unsigned int p = 0;
    
  // Normally each toroidial group can be displayed in the order
  // received. Except for the last plane where it needs to be
  // adjusted by one group. That is if the streamline started in
  // the "correct" place. This is not always the case so it may be
  // necessary to adjust the winding group location by one.
  unsigned int k, offset;
    
  if( p == adjust_plane )
  {
    k = (j-1 + toroidalWindings) % toroidalWindings;
    if( islands )
      offset = 1;
    else
      offset = 0;
  }
  else
  {
    k = j;
    offset = 0;
  }
    
  unsigned int jj = nplanes * toroidalWindings;
    
  if( color == 3 )
    color_value = jj;
    
  // Loop through each point in toroidial group.
  for(unsigned int i=0; i<nnodes; ++i )
  {
    n1 = omesh->add_node(nodes[p][k][(i+offset)%nnodes]);
    ofield->resize_fdata();

    if( color == 2 )
      color_value = (i*toroidalWindings+j)*nplanes + p;
    else if( color == 5 )
      color_value =  i;
        
    ofield->set_value( color_value, n1 );
  }

  // For a surface add in the first point from the adjacent
  // toroidial group. Otherwise for an island add in the
  // first point from the current toroidal group.
  if( !islands )
  {
    if( flip )
      k = (k-skip+toroidalWindings) % toroidalWindings;
    else
      k = (k+skip) % toroidalWindings;
  }

  unsigned int i = nnodes;

  n1 = omesh->add_node(nodes[p][k][0]);
  ofield->resize_fdata();

  if( color == 2 )
    color_value = (i*toroidalWindings+j)*nplanes + p;
  else if( color == 5 )
    color_value =  i;
    
  ofield->set_value( color_value, n1 );
}


} // End namespace FusionPSE
