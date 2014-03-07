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


#ifndef CORE_ALGORITHMS_DATAIO_STREAMDATA_STREAMACQFILE_H
#define CORE_ALGORITHMS_DATAIO_STREAMDATA_STREAMACQFILE_H 1

//! Datatypes that the algorithm uses
#include <Core/Datatypes/Matrix.h>

//! Base class for algorithm
#include <Core/Algorithms/Util/AlgoBase.h>

//! for Windows support
#include <Core/Algorithms/DataIO/share.h>

namespace SCIRunAlgo {

using namespace SCIRun;

class SCISHARE StreamACQFileAlgo : public AlgoBase
{
  public:
    //! Set defaults
    StreamACQFileAlgo()
    {
      add_filename("filename","");
      add_option("method","frame_indices","lead_indices|lead_weights|frame_indices|frame_weights");
      add_handle("gain_table");

      add_int("num_leads",0,Algo::OUT_E);
      add_int("num_frames",0,Algo::OUT_E);
      add_string("time_stamp","",Algo::OUT_E);
      add_string("label","",Algo::OUT_E);
    }
  
    //! Read the header
    bool read_header();
    
    //! Select certain columns or rows
    bool run(MatrixHandle indices,
             MatrixHandle& output);
    
    //! Select a fixed block of data
    bool run(index_type start, 
             index_type end,
             MatrixHandle& output);
};

} // end namespace SCIRunAlgo

#endif

