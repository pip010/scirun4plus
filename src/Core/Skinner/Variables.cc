//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//
//    File   : Variables.cc
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:01:09 2006
#include <Core/Skinner/Variables.h>

#include <Core/XMLUtil/XMLUtil.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Assert.h>

#include <libxml/xmlreader.h>
#include <libxml/catalog.h>
#include <iostream>

#include <string.h>
#include <stdio.h>



namespace SCIRun {
namespace Skinner {


Variables::Variables(const string &id, Variables *parent)
  : variables_(),
    parent_(parent),
    children_(),
    propagate_(),
    alias_()
{
  if (parent) {
    insert("id", parent->get_id()+"/"+id);
    parent->children_.insert(this);
    owned_by_parent_ = true;
  } else if (!id.empty()) {
    insert("id", id);
  }
}


Variables::Variables(const Variables &copy)
  : variables_(copy.variables_),
    parent_(0),
    children_(),
    propagate_(copy.propagate_),
    alias_(copy.alias_),
    owned_by_parent_(false),
    cached_ints_(copy.cached_ints_),
    cached_bools_(copy.cached_bools_),
    cached_doubles_(copy.cached_doubles_),
    cached_strings_(copy.cached_strings_),
    cached_colors_(copy.cached_colors_)
{
}


Variables::~Variables()
{
  children_t::iterator citer = children_.begin();
  while (citer != children_.end()) {
    Variables *child = *(citer++);
    child->parent_ = 0;
    if (child->owned_by_parent_) {
      // Must delete child after incrementing pointer
      // because code below will be invoked by child
      // and remove iterator
      delete child;
    }
  }

  if (parent_) {
    parent_->children_.erase(this);
  }
}


void
Variables::set_parent(Variables *parent)
{
  parent_ = parent;

  // BJW - this didn't use to be here, but I was finding a case where a node's
  //   parent was being deleted by somebody else, and this destructor wanted to
  //   delete the parent again.  So now we keep track from both ends
  parent->children_.insert(this);
}


void
Variables::merge(const Variables *vars)
{
  if (!vars || vars == this) return;
  alias_.insert(vars->alias_.begin(), vars->alias_.end());
  propagate_.insert(vars->propagate_.begin(), vars->propagate_.end());
  name_value_map_t::const_iterator vbeg = vars->variables_.begin();
  name_value_map_t::const_iterator vend = vars->variables_.end();
  for (; vbeg != vend; ++vbeg) {
    //        variables_[vbeg->first] = vbeg->second;
    pair<name_value_map_t::iterator, bool>
      insert = variables_.insert(*vbeg);
    if (!insert.second) {
      insert.first->second = vbeg->second;
    }

    insert.first->second.cache_index_ = -1;
  }
}


string
Variables::type_to_string(Variables::var_type_e type_enum)
{
  switch (type_enum) {
  case UNKNOWN_E: return "unknown";
  case INT_E: return "int";
  case BOOL_E: return "bool";
  case DOUBLE_E: return "double";
  case STRING_E: return "string";
  case COLOR_E: return "color";
  default: throw "Unknown type enum in type_to_string";
  }
  return ""; // shouldnt reach here
}


Variables::var_type_e
Variables::string_to_type(string type_str)
{
  type_str = string_tolower(type_str);
  if (type_str == "unknown") return UNKNOWN_E;
  if (type_str == "int")     return INT_E;
  if (type_str == "bool")    return BOOL_E;
  if (type_str == "float")   return DOUBLE_E;
  if (type_str == "double")  return DOUBLE_E;
  if (type_str == "string")  return STRING_E;
  if (type_str == "color")   return COLOR_E;
  throw "Invalid variable type string: "+type_str;
  return UNKNOWN_E; // shouldn't reach here
}


string &
Variables::find_alias(string &name)
{
  map<string,string>::iterator aiter = alias_.find(name);
  map<string,string>::iterator aend = alias_.end();
  while (aiter != aend) {
    name = aiter->second;
    aiter = alias_.find(name);
  }
  return name;
}


bool
Variables::propagate(const char *name)
{
  return propagate_.find(name) != propagate_.end();
}


pair<Variables *, Variables::value_t *>
Variables::find_value_ptr(const string &inname)
{
  string name = inname;
  Variables *vars_ptr = this;
  Variables::name_value_map_t::iterator iter;

  while (vars_ptr) {
    if (this == vars_ptr || vars_ptr->propagate(name.c_str())) {
      vars_ptr->find_alias(name);
      iter = vars_ptr->variables_.find(name);
      if (iter != vars_ptr->variables_.end()) break;
    }
    vars_ptr = vars_ptr->parent_;
  }

  Variables::value_t *value_ptr = 0;;
  if (vars_ptr) {
    value_ptr = &(iter->second);
  }

  return make_pair(vars_ptr,value_ptr);
}


pair<Variables *, Variables::value_t *>
Variables::insert_variable(const string &name,
                           var_type_e var_type,
                           bool propagate)
{
  // Then it must be inserted
  pair<name_value_map_t::iterator, bool> result =
    variables_.insert(make_pair(name, value_t(name, "", var_type)));
  if (propagate) propagate_.insert(name.c_str());

  return make_pair(this, &(result.first->second));

}


void
Variables::insert(const string &name,
                  const string &in_value,
                  const string &type_str,
                  bool propagate)
{
  string string_value = in_value;
  if (string_value.size() && string_value[0] == '$') {
    string newname = string_value.substr(1,string_value.size()-1);
    alias_[name] = newname;
    if (propagate) propagate_.insert(name.c_str());
    return;
  }

  var_value_t var_val =
    insert_variable(name, string_to_type(type_str), propagate);
  var_val.second->string_value_ = string_value;

  if (!var_val.second->update_cache_from_string(var_val.first))
    throw "invalid conversion";
}


bool
Variables::value_t::update_cache_from_string(Variables *vars)
{
  // Now, Convert the string into the appropriate type,
  // and insert the typed value into the appropriate typed cache vector
  bool success = false;
  switch (var_type_)
  {

  case UNKNOWN_E: {
    success = true; // Nothing to do, its not typed yet
  } break;

  case INT_E: {
    int typed_value;
    if ((success = string_to_int(string_value_,typed_value)))
    {
      vars->set_by_idx(cache_index_, typed_value);
    }
  }  break;

  case BOOL_E: {
    string str = string_toupper(string_value_);
    if (str=="0" || str=="F" || str=="OFF" || str=="FALSE" || str=="NO") {
      vars->set_by_idx(cache_index_, false);
      success = true;
    } else if (str=="1"||str=="T"||str=="ON"||str=="TRUE"||str=="YES") {
      vars->set_by_idx(cache_index_, true);
      success = true;
    }
  } break;

  case DOUBLE_E: {
    double typed_value;
    if ((success = string_to_double(string_value_, typed_value))) {
      vars->set_by_idx(cache_index_, typed_value);
    }
  } break;

  case STRING_E: {
    vars->set_by_idx(cache_index_, string_value_);
    success = true;
  } break;

  case COLOR_E: {
    string str = string_toupper(string_value_);
    if (!str.empty() && str[0] == '#') {
      Color typed_value;
      int r,g,b,a;
      str = str.substr(1, str.length()-1);
      sscanf(str.c_str(), "%02X%02X%02X%02X", &r, &g, &b, &a);

      const double scale = 1.0/255.0;
      typed_value.r = r*scale;
      typed_value.g = g*scale;
      typed_value.b = b*scale;
      typed_value.a = a*scale;

      success = true;
      vars->set_by_idx(cache_index_, typed_value);
    }
  } break;

  default: break;
  }

  cache_current_ = success;
  return success;
}


void
Variables::value_t::update_string_from_cache(Variables *vars)
{
  switch (var_type_)
  {

  case UNKNOWN_E: {
    ASSERT(cache_index_ == -1);
  } break;

  case INT_E: {
    int typed_value;
    vars->get_by_idx(cache_index_, typed_value);
    string_value_ = to_string(typed_value);
  } break;

  case BOOL_E: {
    bool typed_value;
    vars->get_by_idx(cache_index_, typed_value);
    string_value_ = (typed_value ? "1" : "0");
  } break;

  case DOUBLE_E: {
    double typed_value, str_value;
    vars->get_by_idx(cache_index_, typed_value);
    if (string_to_double(string_value_, str_value) &&
        str_value == typed_value) return;
    string_value_ = to_string(typed_value);
  } break;

  case STRING_E: {
    string typed_value;
    vars->get_by_idx(cache_index_, typed_value);
    string_value_ = typed_value;
  } break;

  case COLOR_E: {
    Color typed_value;
    vars->get_by_idx(cache_index_, typed_value);
    char temp[128];
    sprintf(temp, "#%02X%02X%02X%02X",
            (unsigned char)(typed_value.r * 255.0 + 0.5),
            (unsigned char)(typed_value.g * 255.0 + 0.5),
            (unsigned char)(typed_value.b * 255.0 + 0.5),
            (unsigned char)(typed_value.a * 255.0 + 0.5));
    string_value_ = temp;
  } break;

  default: {
    ASSERT(0); // shouldnt reach here
    break;
  }
  }
  cache_current_ = false;
}


int
Variables::get_int(const string &name)
{
  return Var<int>(this,name)();
}


bool
Variables::get_bool(const string &name)
{
  if (!exists(name)) return false;
  Var<bool>val(this,name);
  return val();
}


double
Variables::get_double(const string &name)
{
  return Var<double>(this,name)();
}


string
Variables::get_string(const string &name)
{
  string val;

  pair<Variables *, value_t *> value_ptr = find_value_ptr(name);

  if (!value_ptr.second) {
    throw "get_string failed: "+name;
  }

  if (value_ptr.first &&
      value_ptr.second->var_type_ == STRING_E &&
      value_ptr.second->cache_index_ != -1) {
    val = value_ptr.first->cached_strings_[value_ptr.second->cache_index_];
  } else {
    value_ptr.second->update_string_from_cache(value_ptr.first);
    val = value_ptr.second->string_value_;
  }
  return val;
}


Color
Variables::get_color(const string &name)
{
  return Var<Color>(this,name)();
}


string
Variables::get_id()
{
  return get_string("id");
}


Variables::var_type_e
Variables::get_type_e(const string &name)
{
  var_value_t ptr = find_value_ptr(name);
  if (ptr.second) return ptr.second->var_type_;
  return UNKNOWN_E;
}


Variables::value_t::value_t(string name,
                            string string_value,
                            var_type_e var_type) :
  name_(name),
  string_value_(string_value),
  var_type_(var_type),
  cache_index_(-1)
{
}


bool
Variables::exists(const string &varname)
{
  var_value_t test = find_value_ptr(varname);
  return test.second;
}


void
Variables::set_by_idx(int & idx, const int &value)
{
  set_typed_cache_value<int>(cached_ints_, idx, value);
}


void
Variables::set_by_idx(int & idx, const double &value)
{
  set_typed_cache_value<double>(cached_doubles_, idx, value);
}


void
Variables::set_by_idx(int & idx, const bool &value)
{
  set_typed_cache_value<bool>(cached_bools_, idx, value);
}


void
Variables::set_by_idx(int & idx, const string &value)
{
  set_typed_cache_value<string>(cached_strings_, idx, value);
}


void
Variables::set_by_idx(int & idx, const Color &value)
{
  set_typed_cache_value<Color>(cached_colors_, idx, value);
}


void
Variables::get_by_idx(int & idx, int &value)
{
  ASSERT(idx < (int)cached_ints_.size());
  if (idx < 0) throw "get_by_idx failed on uninitialized value";
  value = cached_ints_[idx];
}


void
Variables::get_by_idx(int & idx, double &value)
{
  ASSERT(idx < (int)cached_doubles_.size());
  if (idx < 0) throw "get_by_idx failed on uninitialized value";
  value = cached_doubles_[idx];
}


void
Variables::get_by_idx(int & idx, bool &value)
{
  ASSERT(idx < (int)cached_bools_.size());
  if (idx < 0) throw "get_by_idx failed on uninitialized value";
  value = cached_bools_[idx];
}


void
Variables::get_by_idx(int & idx, string &value)
{
  ASSERT(idx < (int)cached_strings_.size());
  if (idx < 0) throw "get_by_idx failed on uninitialized value";
  value = cached_strings_[idx];
}


void
Variables::get_by_idx(int & idx, Color &value)
{
  ASSERT(idx < (int)cached_colors_.size());
  if (idx < 0) throw "get_by_idx failed on uninitialized value";
  value = cached_colors_[idx];
}


void
Variables::copy_var(const string &fromid, const string &toid)
{
  var_value_t from = find_value_ptr(fromid);
  var_value_t to = find_value_ptr(toid);
  ASSERT(from.first && from.second && to.first && to.second);
  var_type_e vartype = to.second->var_type_;
  if (vartype == UNKNOWN_E) vartype = from.second->var_type_;
  switch (vartype) {
  case INT_E: { Var<int>(this, toid) = Var<int>(this,fromid)(); } break;
  case BOOL_E: { Var<bool>(this, toid) = Var<bool>(this,fromid)(); } break;
  case DOUBLE_E: {
    Var<double>(this, toid) = Var<double>(this,fromid)();
  } break;
  case STRING_E: {
    Var<string>(this, toid) = Var<string>(this,fromid)();
  } break;
  case COLOR_E: {
    Var<Color>(this, toid) = Var<Color>(this,fromid)();
  } break;
  case UNKNOWN_E: {
    to.second->string_value_ = from.second->string_value_;
  } break;
  default: throw "unknown type in copy_var"; break;
  }
}


bool
Variables::set_by_string(const string &var, const string &val)
{
  if (!exists(var)) return false;
  switch (get_type_e(var)) {
  case INT_E: {
    int temp = 0;
    if (!string_to_int(val, temp)) return false;
    Var<int>(this, var) = temp;
    return true;
  } break;

  case DOUBLE_E: {
    double temp = 0;
    if (!string_to_double(val, temp)) return false;
    pair<Variables *, value_t *> value_ptr = find_value_ptr(var);
    value_ptr.second->string_value_ = val;
    Var<double>(this, var) = temp;
    return true;
  } break;

  case STRING_E: {
    Var<string>(this, var) = val;
    return true;
  } break;

  case BOOL_E: {
    string str = string_toupper(val);
    bool result = false;
    if (str=="0" || str=="F" || str=="OFF" || str=="FALSE" || str=="NO") {
      result = false;
    } else if (str=="1"||str=="T"||str=="ON"||str=="TRUE"||str=="YES") {
      result = true;
    }
    else return false;
    Var<bool>(this, var) = result;
    return true;
  } break;

  // Untested
  case COLOR_E: {
    string str = string_toupper(val);
    if (!str.empty() && str[0] == '#') {
      Color typed_value;
      int r,g,b,a;
      str = str.substr(1, str.length()-1);
      sscanf(str.c_str(), "%02X%02X%02X%02X", &r, &g, &b, &a);

      const double scale = 1.0/255.0;
      typed_value.r = r*scale;
      typed_value.g = g*scale;
      typed_value.b = b*scale;
      typed_value.a = a*scale;

      Var<Color>(this, var) = typed_value;
      return true;
    }
    else return false;
    break;
  }

  case UNKNOWN_E: {
    var_value_t ptr = find_value_ptr(var);
    ASSERT(ptr.second);
    ptr.second->string_value_ = val;
    return true;
  } break;

  default:  // TODO: BOOL_E and COLOR_E
    break;
  }
  return false;
}


}
}
