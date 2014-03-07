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


//    File   : GetFieldZeroCrossings.h
//    Author : Allen Sanderson
//             School of Computing
//             University of Utah
//    Date   : April 2005

#if !defined(GetFieldZeroCrossings_h)
#define GetFieldZeroCrossings_h

#include <Core/Containers/Handle.h>
#include <Core/Datatypes/Field.h>
#include <Core/Util/TypeDescription.h>
#include <Core/Util/DynamicLoader.h>
#include <Core/Math/MiscMath.h>

namespace FusionPSE {

using namespace SCIRun;

class GetFieldZeroCrossingsAlgo : public DynamicAlgoBase
{
public:
  virtual FieldHandle execute(FieldHandle& src, int axis) = 0;
  
  //! support the dynamically compiled algorithm concept
  static CompileInfoHandle get_compile_info(const TypeDescription *iftd,
					    const string &sftd);
};

template< class IFIELD, class OFIELD >
class GetFieldZeroCrossingsAlgoT : public GetFieldZeroCrossingsAlgo
{
public:
  //! virtual interface. 
  virtual FieldHandle execute(FieldHandle& src, int axis);
};


template< class IFIELD, class OFIELD >
FieldHandle
GetFieldZeroCrossingsAlgoT<IFIELD, OFIELD>::execute(FieldHandle& ifield_h,
						    int axis)
{
  IFIELD *ifield = (IFIELD *) ifield_h.get_rep();
  OFIELD *ofield = new OFIELD(ifield->get_typed_mesh());

  typename OFIELD::mesh_handle_type omesh = ofield->get_typed_mesh();

  typename OFIELD::fdata_type::iterator in  = ofield->fdata().begin();
  typename OFIELD::fdata_type::iterator end = ofield->fdata().end();

  // Zero out the output data.
  while (in != end)
  {
    *in = 0.0;
    ++in;
  }

  // Nodes are working. Some mesh types need this.
  omesh->synchronize(Mesh::NODES_E);
  if (omesh->dimensionality() >= 1) omesh->synchronize(Mesh::EDGES_E);
  if (omesh->dimensionality() >= 2) omesh->synchronize(Mesh::FACES_E);
  if (omesh->dimensionality() >= 3) omesh->synchronize(Mesh::CELLS_E);

  typename OFIELD::mesh_type::Node::array_type nodes;

  typename OFIELD::mesh_type::Edge::iterator iEdgeItr;
  typename OFIELD::mesh_type::Edge::iterator oEdgeItr;

  omesh->begin( iEdgeItr );
  omesh->end  ( oEdgeItr );

  Point p0, p1;

  // Iterate through each edge.
  while( iEdgeItr != oEdgeItr ) {

    omesh->get_nodes(nodes, *iEdgeItr);

    // Get the point and value at the locations.
    omesh->get_center(p0, nodes[0]);
    omesh->get_center(p1, nodes[1]);

    for( unsigned int i=0; i<2; i++ ) {

      Vector diff;

      // Get the derivative between the two points.
      if( i == 0 ) diff = (Vector) p0 - (Vector) p1;
      else         diff = (Vector) p1 - (Vector) p0;

      // Normalize to be -1 or +1
      if( diff.x() > 0 ) diff.x( 1 ); else if( diff.x() < 0 ) diff.x( -1 );
      if( diff.y() > 0 ) diff.y( 1 ); else if( diff.y() < 0 ) diff.y( -1 );
      if( diff.z() > 0 ) diff.z( 1 ); else if( diff.z() < 0 ) diff.z( -1 );

      // Sum the derivative sign at each point. If the sum is zero there
      // is a zero crossing.

      // Get the current sum.
      double val = ofield->value( nodes[i] );
      
      // Add the new edge to the sum.
      if( axis == 0 )      val += diff.x();
      else if( axis == 1 ) val += diff.y();
      else if( axis == 2 ) val += diff.z();

      // Resave the sum.
      ofield->set_value(val, nodes[i]);
    }

    ++iEdgeItr;
  }

  // For a zero crossing to occur at least two edges must go in the
  // same direction. So those that are -2 or 2 become 1 (true)
  // otherwise it becomes 0 (false).

  in  = ofield->fdata().begin();
  end = ofield->fdata().end();

  while (in != end)
  {
    if( fabs( *in ) > 1.0 )
      *in = 1.00001; // Make this just over 1.0 because of rounding.
    else
      *in = 0.0;

    ++in;
  }


  typename IFIELD::fdata_type::iterator dataItr = ifield->fdata().begin();
  typename IFIELD::mesh_type::Node::iterator inodeItr, inodeEnd;

  omesh->begin( inodeItr );
  omesh->end( inodeEnd );
  dataItr = ofield->fdata().begin();

  Point currPt;

  omesh->get_center(currPt, *inodeItr);

  Point lastPt = currPt;
  Point penultimatePt = lastPt;

  bool check = true;

  while (inodeItr != inodeEnd) {
    // Get the current point.
    omesh->get_center(currPt, *inodeItr);

      // Find the positive zero crossing.
    if( check &&
	0.0 < currPt.z() &&
	penultimatePt.z() <= lastPt.z() &&
	currPt.z() < lastPt.z() ) {
      *dataItr = 1.001;

      check = false;
    } else {
      *dataItr = 0;
    }

    if( currPt.z() < 0.0 && check == false )
      check = true;
    
    penultimatePt = lastPt;
    lastPt  = currPt;

    ++inodeItr;
    ++dataItr;
  }


  return FieldHandle(ofield);
}

} // end namespace FusionPSE

#endif // GetFieldZeroCrossings_h
