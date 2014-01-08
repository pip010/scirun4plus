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

#ifndef CORE_ALGORITHMS_UTIL_ALGOCALLBACK_H
#define CORE_ALGORITHMS_UTIL_ALGOCALLBACK_H 1

#include <list>

//! For Windows
#include <Core/Algorithms/Util/share.h>

namespace SCIRunAlgo {

//! Base class for writing callbacks. Group data and the callback function
//! into a class based on this one class.

class SCISHARE AlgoCallBack {

  public:
    virtual bool callback() = 0;
    
    virtual ~AlgoCallBack() {}
};


//! Class that the algorithm derives from

class SCISHARE AlgoCallBackList {

  public:
  
    //! Add a callback to the stack
    void add_callback(AlgoCallBack* callback);
      
    //! Remove a callback from the stack  
    void remove_callback(AlgoCallBack* callback);
    
  //! These functions are called by the algorithm
  public:
  
    //! Execute all the callbacks
    void do_callbacks();
  
    //! Are there any callback.
    //! Mainly to establish whether data needs to be synchronized before
    //! doing callbacks.
    bool have_callbacks();

  private:
    
    std::list<AlgoCallBack*>   callbacks_;

};

} //namespace SCIRun

#endif
