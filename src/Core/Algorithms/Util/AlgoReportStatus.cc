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

#include <Core/Algorithms/Util/AlgoReportStatus.h>

namespace SCIRunAlgo {


AlgoReportStatus::AlgoReportStatus()
{ 
  pr_ = &(default_pr_); 
}

void 
AlgoReportStatus::set_progress_reporter(SCIRun::ProgressReporter* pr)
{ 
  pr_ = pr;  
}

void 
AlgoReportStatus::set_progress_reporter(AlgoReportStatus* algo)
{ 
  pr_ = algo->get_progress_reporter();  
}
  
SCIRun::ProgressReporter* 
AlgoReportStatus::get_progress_reporter() const
{ 
  return (pr_); 
}
    

void 
AlgoReportStatus::error(const std::string& error) const
{ 
  pr_->error(error); 
}

void 
AlgoReportStatus::warning(const std::string& warning) const
{ 
  pr_->warning(warning); 
}

void 
AlgoReportStatus::remark(const std::string& remark) const
{ 
  pr_->remark(remark);
}

void 
AlgoReportStatus::status(const std::string& status) const
{ 
  pr_->status(status);
}

void 
AlgoReportStatus::update_progress(double percent) const
{ 
  pr_->update_progress(percent); 
}

void 
AlgoReportStatus::update_progress(SCIRun::index_type current, SCIRun::size_type max) const
{ 
  pr_->update_progress(current,max); 
}

void 
AlgoReportStatus::algo_start(const std::string& tag)
{
  pr_->report_start(tag);
}

void
AlgoReportStatus::algo_end()
{
  pr_->report_end();
}

bool
AlgoReportStatus::check_abort() const
{
  // TODO : Add this hook to progress reporter
  return (false);
}

ScopedAlgorithmReporter::ScopedAlgorithmReporter(AlgoReportStatus* algo, const std::string& tag) : algo_(algo)
{
  algo_->algo_start(tag);
}

ScopedAlgorithmReporter::~ScopedAlgorithmReporter()
{
  algo_->algo_end();
}

} //namespace SCIRunAlgo
