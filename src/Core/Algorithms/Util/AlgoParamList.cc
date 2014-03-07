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

#include <boost/algorithm/string.hpp>
#include <Core/Algorithms/Util/AlgoParamList.h>
#include <Core/Util/MemoryUtil.h>
#include <ctype.h>

namespace SCIRunAlgo {

bool 
AlgoParamList::set_bool(const std::string& key, bool value)
{
  BoolParam* param = dynamic_cast<BoolParam*>(parameters_[key]);
  if (param) { param->value() = value; return (true); }
  
  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}

bool
AlgoParamList::get_bool(const std::string& key, bool& value)
{
  BoolParam* param = dynamic_cast<BoolParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}


bool
AlgoParamList::get_bool(const std::string& key)
{
  BoolParam* param = dynamic_cast<BoolParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

bool 
AlgoParamList::set_int(const std::string& key, int value)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) 
  { 
    if ((!param->limit()) ||
        (param->limit() && ((value < param->min())||(value > param->max()))))
    {
      param->value() = value; return (true); 
    }
    throw std::string("parameter \""+key+"\" is out of limits");
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}

bool
AlgoParamList::get_int(const std::string& key, int& value)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

int
AlgoParamList::get_int(const std::string& key)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (0); 
  }
}


bool
AlgoParamList::get_int_limits(const std::string& key, int& min, int& max)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) { min = param->min(); max = param->max(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}


bool 
AlgoParamList::set_index(const std::string& key, index_type value)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) 
  { 
    if ((!param->limit()) ||
        (param->limit() && ((value < param->min())||(value > param->max()))))
    {
      param->value() = value; return (true); 
    }
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}

bool
AlgoParamList::get_index(const std::string& key, index_type& value)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

index_type
AlgoParamList::get_index(const std::string& key)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (0); 
  }
}


bool
AlgoParamList::get_index_limits(const std::string& key, index_type& min, index_type& max)
{
  IntegerParam* param = dynamic_cast<IntegerParam*>(parameters_[key]);
  if (param) { min = param->min(); max = param->max(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}




bool 
AlgoParamList::set_scalar(const std::string& key, double value)
{
  ScalarParam* param = dynamic_cast<ScalarParam*>(parameters_[key]);
  if (param) 
  { 
    if ((!param->limit()) ||
        (param->limit() && ((value < param->min())||(value > param->max()))))
    {
      param->value() = value; return (true); 
    }
    throw std::string("parameter \""+key+"\" is out of limits");
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}

bool
AlgoParamList::get_scalar(const std::string& key, double& value)
{
  ScalarParam* param = dynamic_cast<ScalarParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

double
AlgoParamList::get_scalar(const std::string& key)
{
  ScalarParam* param = dynamic_cast<ScalarParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (0.0); 
  }
}

bool
AlgoParamList::get_scalar_limits(const std::string& key, double& min, double& max)
{
  ScalarParam* param = dynamic_cast<ScalarParam*>(parameters_[key]);
  if (param) { min = param->min(); max = param->max(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

bool
AlgoParamList::set_color(const std::string& key, SCIRun::Color value)
{
  ColorParam* param = dynamic_cast<ColorParam*>(parameters_[key]);
  if (param) 
  {
    param->value() = value;
    return (true); 
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}


bool
AlgoParamList::get_color(const std::string& key, SCIRun::Color& value)
{
  ColorParam* param = dynamic_cast<ColorParam*>(parameters_[key]);
  if (param) {
    value = param->value();
    return (true);
  }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

SCIRun::Color
AlgoParamList::get_color(const std::string& key)
{
  ColorParam* param = dynamic_cast<ColorParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  {
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (SCIRun::Color(0.0,0.0,0.0)); 
  }
}

bool
AlgoParamList::set_colormap(const std::string& key, SCIRun::ColorMap* value)
{
  ColorMapParam* param = dynamic_cast<ColorMapParam*>(parameters_[key]);
  if (param)
  {
    param->value() = value;
    return (true); 
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}

bool
AlgoParamList::get_colormap(const std::string& key, SCIRun::ColorMap*& value)
{
  ColorMapParam* param = dynamic_cast<ColorMapParam*>(parameters_[key]);
  if (param)
  {
    value = param->value().get_rep();
    return (true);
  }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

SCIRun::ColorMap*
AlgoParamList::get_colormap(const std::string& key)
{
  ColorMapParam* param = dynamic_cast<ColorMapParam*>(parameters_[key]);
  if (param) { return(param->value().get_rep()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (0); 
  }
}

bool
AlgoParamList::set_vector(const std::string& key, SCIRun::Vector value)
{
  VectorParam* param = dynamic_cast<VectorParam*>(parameters_[key]);
  if (param) 
  {
    param->value() = value; return (true); 
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}

bool
AlgoParamList::get_vector(const std::string& key, SCIRun::Vector& value)
{
  VectorParam* param = dynamic_cast<VectorParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

SCIRun::Vector
AlgoParamList::get_vector(const std::string& key)
{
  VectorParam* param = dynamic_cast<VectorParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (SCIRun::Vector(0.0,0.0,0.0)); 
  }
}


bool
AlgoParamList::set_point(const std::string& key, SCIRun::Point value)
{
  PointParam* param = dynamic_cast<PointParam*>(parameters_[key]);
  if (param) 
  {
    param->value() = value; return (true); 
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false); 
}

bool
AlgoParamList::get_point(const std::string& key, SCIRun::Point& value)
{
  PointParam* param = dynamic_cast<PointParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

SCIRun::Point
AlgoParamList::get_point(const std::string& key)
{
  PointParam* param = dynamic_cast<PointParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (SCIRun::Point(0.0,0.0,0.0)); 
  }
}

bool
AlgoParamList::set_option(const std::string& key, const std::string& value)
{
  OptionParam* param = dynamic_cast<OptionParam*>(parameters_[key]);
  std::string valueLower = boost::to_lower_copy(value);
  for (size_t j=0; j<valueLower.size(); j++) valueLower[j] = tolower(valueLower[j]);
  
  if (param) 
  {
    std::vector<std::string>::iterator it = param->options().begin();
    std::vector<std::string>::iterator it_end = param->options().end();
    while (it != it_end)
    {
      if ((*it) == valueLower) { param->value() = valueLower; return (true); }
      ++it;
    }
    throw std::string("parameter \""+key+"\" has no option \""+valueLower+"\"");    
  }

  throw std::string("key \""+key+"\" was not defined in algorithm");
  return (false);
}

bool
AlgoParamList::get_option(const std::string& key, std::string& value)
{
  OptionParam* param = dynamic_cast<OptionParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
   throw std::string("key \""+key+"\" was not defined in algorithm");
   return (false); 
  }
}

std::string
AlgoParamList::get_option(const std::string& key)
{
  OptionParam* param = dynamic_cast<OptionParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (""); 
  }
}

bool
AlgoParamList::check_option(const std::string& key,const std::string& value)
{
  std::string valueLower = boost::to_lower_copy(value);
  OptionParam* param = dynamic_cast<OptionParam*>(parameters_[key]);
  if (param) { return(param->value() == valueLower); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }  
}


bool 
AlgoParamList::set_string(const std::string& key, const std::string& value)
{
  StringParam* param = dynamic_cast<StringParam*>(parameters_[key]);
  if (param) { param->value() = value; return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

bool
AlgoParamList::get_string(const std::string& key, std::string& value)
{
  StringParam* param = dynamic_cast<StringParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

std::string
AlgoParamList::get_string(const std::string& key)
{
  StringParam* param = dynamic_cast<StringParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (""); 
  }
}

bool 
AlgoParamList::set_filename(const std::string& key, const std::string& value)
{
  FilenameParam* param = dynamic_cast<FilenameParam*>(parameters_[key]);
  if (param) { param->value() = value; return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

bool
AlgoParamList::get_filename(const std::string& key, std::string& value)
{
  FilenameParam* param = dynamic_cast<FilenameParam*>(parameters_[key]);
  if (param) { value = param->value(); return (true); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (false); 
  }
}

std::string
AlgoParamList::get_filename(const std::string& key)
{
  FilenameParam* param = dynamic_cast<FilenameParam*>(parameters_[key]);
  if (param) { return(param->value()); }
  else 
  { 
    throw std::string("key \""+key+"\" was not defined in algorithm");
    return (""); 
  }
}


AlgoParamList::~AlgoParamList()
{
  SCIRun::delete_all_values(parameters_);
}

void
AlgoParamList::add_bool(const std::string& key, bool defval, int type)
{
  parameters_[key] = new BoolParam(defval,type);
}

void
AlgoParamList::add_int(const std::string& key, int defval, int type)
{
  parameters_[key] = new IntegerParam(defval,type);
}

void
AlgoParamList::add_int(const std::string& key, int defval, int min, int max, int type)
{
  parameters_[key] = new IntegerParam(defval,min,max,type);
}

void
AlgoParamList::add_index(const std::string& key, index_type defval, int type)
{
  parameters_[key] = new IntegerParam(defval,type);
}

void
AlgoParamList::add_index(const std::string& key, index_type defval, index_type min, index_type max, int type)
{
  parameters_[key] = new IntegerParam(defval,min,max,type);
}

void
AlgoParamList::add_scalar(const std::string& key, double defval, int type)
{
  parameters_[key] = new ScalarParam(defval,type);
}


void
AlgoParamList::add_color(const std::string& key, SCIRun::Color defval, int type)
{
  parameters_[key] = new ColorParam(defval,type);
}

void
AlgoParamList::add_vector(const std::string& key, SCIRun::Vector defval, int type)
{
  parameters_[key] = new VectorParam(defval,type);
}

void
AlgoParamList::add_point(const std::string& key, SCIRun::Point defval, int type)
{
  parameters_[key] = new PointParam(defval,type);
}

void
AlgoParamList::add_colormap(const std::string& key, SCIRun::ColorMap* defval, int type)
{
  parameters_[key] = new ColorMapParam(defval,type);
}

void
AlgoParamList::add_scalar(const std::string& key, double defval, double min, double max, int type)
{
  parameters_[key] = new ScalarParam(defval,min,max,type);
}

void
AlgoParamList::add_option(const std::string& key, const std::string& defval, const std::string& options, int type)
{
  std::vector<std::string> opts;
  std::string optionsLower = boost::to_lower_copy(options);
  while(true)
  {
    size_t loc = optionsLower.find('|');
    if (loc >= optionsLower.size())
    {
      opts.push_back(optionsLower);
      break;
    }
    opts.push_back(optionsLower.substr(0,loc));
    optionsLower = optionsLower.substr(loc+1);
  }
  
  parameters_[key] = new OptionParam(defval,opts,type);
}

void
AlgoParamList::add_string(const std::string& key, const std::string& defval, int type)
{
  parameters_[key] = new StringParam(defval,type);
}

void
AlgoParamList::add_filename(const std::string& key, const std::string& defval, int type)
{
  parameters_[key] = new FilenameParam(defval,type);
}

void
AlgoParamList::add_handle(const std::string& key, int type)
{
  parameters_[key] = new HandleParam(type);
}

}
