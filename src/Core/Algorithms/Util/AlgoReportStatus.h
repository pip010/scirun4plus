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

#ifndef CORE_ALGORITHMS_UTIL_ALGOREPORTSTATUS_H
#define CORE_ALGORITHMS_UTIL_ALGOREPORTSTATUS_H 1

#include <Core/Datatypes/Types.h>
#include <Core/Util/ProgressReporter.h>

//! For Windows
#include <Core/Algorithms/Util/share.h>

namespace SCIRunAlgo {

class AlgoReportStatus;

class SCISHARE AlgoReportStatus {

  public:
  
    AlgoReportStatus();
  
    //! Functions for setting progress reporter
    //! The progress reporter relays information to other parts of the program.
    void set_progress_reporter(SCIRun::ProgressReporter* pr);
    void set_progress_reporter(AlgoReportStatus* algo);
    SCIRun::ProgressReporter* get_progress_reporter() const;
    
    //! functions for the algorithm, so it can forward errors if needed
    void error(const std::string& error) const;
    void warning(const std::string& warning) const;
    void remark(const std::string& remark) const;
    void status(const std::string& status) const;

    //! Tag algorithms so we can relay this information to the user
    //! These functions tag which algorithm is currently running
    void algo_start(const std::string& tag);
    void algo_end();
    
    //! Report progress of the algorithm
    void update_progress(double percent) const;
    void update_progress(SCIRun::index_type current, SCIRun::size_type max) const;

    //! Check abort status, by polling this function
    bool check_abort() const;
    
  private:
    // ProgressReporter
    SCIRun::ProgressReporter*   pr_;
    SCIRun::ProgressReporter    default_pr_;
    
};

class SCISHARE ScopedAlgorithmReporter
{
public:
  ScopedAlgorithmReporter(AlgoReportStatus* algo, const std::string& tag);
  ~ScopedAlgorithmReporter();
private:
  AlgoReportStatus* algo_;
};

} //namespace SCIRunAlgo

#endif
