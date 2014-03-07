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
 *   September 2005
 *
 */

#if !defined(StreamlineAnalyzer_h)
#define StreamlineAnalyzer_h

#include <Packages/FusionPSE/Dataflow/Modules/Fields/StreamlineAnalyzerLib.h>
#include <Core/Datatypes/Field.h>
#include <sstream>

namespace FusionPSE {

using namespace SCIRun;

class StreamlineAnalyzerAlgo : public DynamicAlgoBase
{
  public:

    //! virtual interface. 
    virtual void execute(FieldHandle& slsrc,
             FieldHandle& dst,
             FieldHandle& pccsrc,
             FieldHandle& pccdst,
             FieldHandle& pcssrc,
             FieldHandle& pcsdst,
             vector< double > &planes,
             unsigned int color,
             unsigned int showIslands,
             unsigned int islandCentroids,
             unsigned int overlaps,
             unsigned int maxToroidalWinding,
             unsigned int override,
             float hitrate,
             vector< pair< unsigned int,
             unsigned int > > &topology ) = 0;
    
  protected:
    FieldlineLib FLlib;
};


} // End namespace FusionPSE

#endif // StreamlineAnalyzer_h
