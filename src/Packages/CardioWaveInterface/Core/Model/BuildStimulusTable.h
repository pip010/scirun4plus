/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
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

#ifndef CARDIOWAVE_CORE_FIELDS_BUILDSTIMULUSTABLE_H
#define CARDIOWAVE_CORE_FIELDS_BUILDSTIMULUSTABLE_H 1

#include <Core/Util/ProgressReporter.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Matrix.h>

#include <vector>

namespace CardioWaveInterface {

using namespace SCIRun;

class stimulusparam_type {
  public:
    Field::index_type node;
    double            weight;
};


typedef std::vector<stimulusparam_type> StimulusTable;
typedef std::vector<stimulusparam_type> ReferenceTable;

inline bool operator==(const stimulusparam_type& p1,const stimulusparam_type& p2)
{
  if (p1.node == p2.node) return (true);
  return (false);
}    

inline bool operator<(const stimulusparam_type& p1, const stimulusparam_type& p2)
{
  if (p1.node < p2.node) return(true); 
  return (false);
}


class BuildStimulusTableAlgo
{
  public:
    bool BuildStimulusTable(ProgressReporter *pr, FieldHandle ElementType, 
                            FieldHandle Stimulus, MatrixHandle CompToGeom, 
                            double domaintypemin, double domaintypemax, 
                            bool selectbynode, StimulusTable& stimulustablelist);
};

} // end namespace ModelCreation

#endif 

