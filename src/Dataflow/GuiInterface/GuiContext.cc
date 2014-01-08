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

/*
 *  GuiContext.cc:
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 2002
 *
 */

#include <Core/Util/StringUtil.h>
#include <Core/Util/Environment.h>

#include <Core/Math/MiscMath.h>
#include <Core/Util/Assert.h>

#include <Dataflow/GuiInterface/GuiContext.h>
#include <Dataflow/GuiInterface/TCLInterface.h>

#include <iostream>
#include <sstream>
#include <iomanip>

using namespace SCIRun;

GuiContext::GuiContext(const std::string& name,
                       bool save, // default = true
                       GuiContext *parent) // default = 0
  : parent_(parent),
    name_(name),
    children_(),
    context_state_(CACHE_E),
    active_(true)
{
  context_state_ |= save?SAVE_E:0;
  if (save) tcl_setVarStates();
}


GuiContext::~GuiContext()
{  
  if (parent_) 
  {
    for(std::vector<GuiContext*>::iterator iter = parent_->children_.begin();
          iter != parent_->children_.end(); ++iter) 
    {
      if (*iter == this)
      {
        parent_->children_.erase(iter);
        return;
      } 
    }
  }
}



GuiContext* 
GuiContext::subVar(const std::string& subname, bool saveChild)
{
  dontSave(); // Do not save intermediate nodes

  GuiContext* child = find_child(name_+"-"+subname);
  if (child == 0)
  {
    child = new GuiContext(name_+"-"+subname, saveChild, this);
    children_.push_back(child);
  }
  return child;
}


GuiContext * 
GuiContext::find_child(const std::string &subname)
{
  std::string fullname = name_+"-"+subname;
  for (unsigned int c = 0; c < children_.size(); ++c)
  {
    if (fullname == children_[c]->name_)
    {
      return children_[c];
    }
  }
  return 0;
}

void
GuiContext::inactivate()
{
  active_ = false;
  for (unsigned int c = 0; c < children_.size(); ++c) children_[c]->inactivate();
}


bool
GuiContext::get(std::string& value)
{
  if (context_state_ & CACHED_E) return true;
  if(!TCLInterface::get(name_, value,this)) return false;

  context_state_ |= CACHED_E;
  return true;
}

void
GuiContext::async_get(std::string* value)
{
  if (context_state_ & CACHED_E) return;
  TCLInterface::async_get(name_, value,this);
  context_state_ |= CACHED_E;
}

void
GuiContext::set(const std::string& value)
{
  context_state_ &= ~CACHED_E;
  TCLInterface::set(name_, value,this);
}

void
GuiContext::async_set(const std::string& value)
{
  context_state_ &= ~CACHED_E;
  TCLInterface::async_set(name_, value,this);
}


bool
GuiContext::get(double& value)
{
  if (context_state_ & CACHED_E) return true;
  if(!TCLInterface::get(name_, value,this)) return false;

  context_state_ |= CACHED_E;
  return true;
}

void
GuiContext::async_get(double* value)
{
  if (context_state_ & CACHED_E) return;
  TCLInterface::async_get(name_, value,this);
  context_state_ |= CACHED_E;
}

void
GuiContext::set(double value)
{
  context_state_ &= ~CACHED_E;
  TCLInterface::set(name_, value, this);
}

void
GuiContext::async_set(double value)
{
  context_state_ &= ~CACHED_E;
  TCLInterface::async_set(name_, value, this);
}


bool
GuiContext::get(int& value)
{
  if (context_state_ & CACHED_E) return true;
  if(!TCLInterface::get(name_, value,this)) return false;

  context_state_ |= CACHED_E;
  return true;
}

void
GuiContext::async_get(int* value)
{
  if (context_state_ & CACHED_E) return;
  TCLInterface::async_get(name_, value,this);
  context_state_ |= CACHED_E;
}

void
GuiContext::set(int value)
{
  context_state_ &= ~CACHED_E;
  TCLInterface::set(name_, value,this);
}


void
GuiContext::async_set(int value)
{
  context_state_ &= ~CACHED_E;
  TCLInterface::async_set(name_, value,this);
}


void
GuiContext::synchronize()
{ 
  TCLInterface::synchronize(); 
}

void
GuiContext::reset()
{
  TCLInterface::lock();
  if (active_ == false)
  {
    TCLInterface::unlock();
    return;
  }
  TCLInterface::unlock();
  
  context_state_ &= ~CACHED_E;
  for(std::vector<GuiContext*>::iterator iter = children_.begin();
        iter != children_.end(); ++iter) (*iter)->reset();
}

std::string
GuiContext::getfullname()
{
  return name_;
}

void
GuiContext::tcl_setVarStates()
{
  const std::string save_flag = (context_state_ & SAVE_E)?"1":"0";
  const std::string sub_flag = (context_state_ & SUBSTITUTE_DATADIR_E)?" 1":" 0";
  const std::string filename_flag = (context_state_ & IS_FILENAME_E)?" 1":" 0";
  TCLInterface::async_execute("setVarStates \""+name_+"\" "+save_flag+sub_flag+filename_flag,this);
}

void
GuiContext::dontSave()
{
  if ((context_state_ & SAVE_E) == 0) return;
  context_state_ &= ~SAVE_E;
  tcl_setVarStates();
}

void
GuiContext::doSave()
{
  if ((context_state_ & SAVE_E) == SAVE_E) return;
  context_state_ |= SAVE_E;
  tcl_setVarStates();
}


void
GuiContext::dontReplaceEnv()
{
  if ((context_state_ & SUBSTITUTE_ENVIRONMENT_E) == 0) return;
  context_state_ &= ~SUBSTITUTE_ENVIRONMENT_E;
}

void
GuiContext::doReplaceEnv()
{
  if ((context_state_ & SUBSTITUTE_ENVIRONMENT_E) == SUBSTITUTE_ENVIRONMENT_E) return;
  context_state_ |= SUBSTITUTE_ENVIRONMENT_E;
}

bool
GuiContext::replace_env()
{
  return (context_state_ & SUBSTITUTE_ENVIRONMENT_E);
}


void
GuiContext::dontSubstituteDatadir()
{
  if ((context_state_ & SUBSTITUTE_DATADIR_E) == 0) return;
  context_state_ &= ~SUBSTITUTE_DATADIR_E;
  tcl_setVarStates();
}

void
GuiContext::doSubstituteDatadir()
{
  if ((context_state_ & SUBSTITUTE_DATADIR_E) == SUBSTITUTE_DATADIR_E) return;
  context_state_ |= SUBSTITUTE_DATADIR_E;
  tcl_setVarStates();
}


void
GuiContext::unsetIsFilename()
{
  if ((context_state_ & IS_FILENAME_E) == 0) return;
  context_state_ &= ~IS_FILENAME_E;
  tcl_setVarStates();
}

void
GuiContext::setIsFilename()
{
  if ((context_state_ & IS_FILENAME_E) == IS_FILENAME_E) return;
  context_state_ |= IS_FILENAME_E;
  tcl_setVarStates();
}


