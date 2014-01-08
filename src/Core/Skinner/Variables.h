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
//    File   : Variables.h
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:03:25 2006

#ifndef Skinner_Variables_H
#define Skinner_Variables_H

#include <Core/Skinner/Color.h>
#include <Core/Util/Assert.h>

#include <string>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <functional>

#include <Core/Skinner/share.h>

namespace SCIRun {
namespace Skinner {

class Variables;

template <class T>
class Var {
private:
  friend class Variables;
  Variables *       scope_;
  int *             scope_index_;

public:
  Var(Variables *vars, const std::string &name);
  Var(Variables *vars, const std::string &name, const T &);
  Var() : scope_(0), scope_index_(0) {}
  Var(const Var<T> &copy) :
    scope_(copy.scope_),
    scope_index_(copy.scope_index_)
  {
  }
  bool exists() { return scope_index_ && (*scope_index_ >= 0); }
  operator T();
  Var<T> & operator= (const T& rhs);
  Var<T> & operator= (const Var<T>& rhs);
  Var<T> & operator|=(const T& rhs);
  Var<T> & operator|=(const Var<T>& rhs);
  T operator()();
};


class SCISHARE Variables
{
public:

  Variables         (const std::string &id="", Variables *parent=0);
  Variables         (const Variables &copy);

  virtual           ~Variables();
  void              merge(const Variables *);
  void              insert(const std::string &name,
                           const std::string &value,
                           const std:: &type_str = "string",
                           bool propagate = false);
  void              copy_var(const std::string &from, const std::string &to);

  bool              exists(const std::string &varname);

  std::string       get_id();
  int               get_int(const std::string &);
  double            get_double(const std::string &);
  bool              get_bool(const std::string &);
  Color             get_color(const std::string &);
  std::string       get_string(const std::string &);
  void              set_parent(Variables *);

  enum var_type_e {
    UNKNOWN_E,
    INT_E,
    BOOL_E,
    DOUBLE_E,
    STRING_E,
    COLOR_E
  };
  var_type_e        get_type_e(const std::string &);
  bool              set_by_string(const std::string &, const std::string &);

  struct SCISHARE value_t {
    value_t(std::string, std::string, var_type_e);
    bool              update_cache_from_string(Variables *);
    void              update_string_from_cache(Variables *);
    std::string            name_;
    std::string            string_value_;
    var_type_e        var_type_;
    int               cache_index_;
    bool              cache_current_;
  };

private:
  // let Var access private data
  template <class T>
  friend class Var;


  struct ltstr
  {
    bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) < 0;
    }
  };

  typedef std::map<std::string, value_t>              name_value_map_t;
  typedef set<Variables *>                  children_t;
  typedef std::pair<Variables *, value_t *>      var_value_t;
  typedef std::map<std::string, std::string>               alias_t;
  typedef set<std::string>          propagate_t;

  var_value_t               insert_variable(const std::string &,
                                            var_type_e,
                                            bool);

  var_value_t               find_value_ptr(const std::string &);


  template<class T>
  int                       set_typed_cache_value(std::vector<T> &, int &,
                                                  const T &);
  std::string &                  find_alias(std::string &);
  bool                      propagate(const char *);
  static var_type_e         string_to_type(std::string);
  static std::string             type_to_string(var_type_e);
  var_type_e                type_to_enum(int) { return INT_E; }
  var_type_e                type_to_enum(bool){ return BOOL_E; }
  var_type_e                type_to_enum(double){ return DOUBLE_E; }
  var_type_e                type_to_enum(std::string){ return STRING_E; }
  var_type_e                type_to_enum(Color){ return COLOR_E; }

  void                      set_by_idx(int &, const int &);
  void                      set_by_idx(int &, const bool &);
  void                      set_by_idx(int &, const double &);
  void                      set_by_idx(int &, const std::string &);
  void                      set_by_idx(int &, const Skinner::Color &);

  void                      get_by_idx(int &, int &);
  void                      get_by_idx(int &, bool &);
  void                      get_by_idx(int &, double &);
  void                      get_by_idx(int &, std::string &);
  void                      get_by_idx(int &, Skinner::Color &);

  name_value_map_t          variables_;
  Variables *               parent_;
  children_t                children_;
  propagate_t               propagate_;
  alias_t                   alias_;

  // we must use this, since Vars can either be owned by other Variables
  // or by something else.
  bool                      owned_by_parent_;

  std::vector<int>               cached_ints_;
  std::vector<bool>              cached_bools_;
  std::vector<double>            cached_doubles_;
  std::vector<std::string>       cached_strings_;
  std::vector<Skinner::Color>    cached_colors_;
};


template <class T>
Var<T>::operator T()
{
  T temp;
  ASSERT(this->scope_index_  && (*this->scope_index_ >= 0));
  this->scope_->get_by_idx(*this->scope_index_, temp);
  return temp;
}


template <class T>
Var<T> &
Var<T>::operator= (const T& rhs)
{
  ASSERT(this->scope_);
  ASSERT(this->scope_index_ && *this->scope_index_ >= -1);
  this->scope_->set_by_idx(*this->scope_index_, rhs);
  return *this;
}


template <class T>
Var<T> &
Var<T>::operator= (const Var<T>& rhs)
{
  this->scope_ = rhs.scope_;
  this->scope_index_ = rhs.scope_index_;
  return *this;
}


template <class T>
Var<T> &
Var<T>::operator|= (const Var<T>& rhs)
{
  if (this->scope_index_ &&
      (*this->scope_index_ == -1) &&
      (*rhs.scope_index_ != -1)) {
    return operator=(rhs);
  }
  return *this;
}


template <class T>
Var<T> &
Var<T>::operator|= (const T& rhs)
{
  ASSERT(this->scope_);
  if (this->scope_index_ && (*this->scope_index_ == -1)) {
    this->scope_->set_by_idx(*this->scope_index_, rhs);
  }
  ASSERT(this->scope_index_  && (*this->scope_index_ >= 0));
  return *this;
}


template <class T>
T
Var<T>::operator()()
{
  T temp;
  ASSERT(this->scope_index_  && (*this->scope_index_ >= 0));
  this->scope_->get_by_idx(*this->scope_index_, temp);
  return temp;
}


template<class T>
int
Variables::set_typed_cache_value(std::vector<T> &cache_vector,
                                 int &index,
                                 const T &typed_value)
{
  int size = int(cache_vector.size());
  ASSERT(index < size);
  if (index >= 0) {
    cache_vector[index] = typed_value;
  } else {
    index = size;
    cache_vector.push_back(typed_value);
  }
  return index;
}


template <class T>
Var<T>::Var(Variables *vars, const std::string &name, const T &init)
{
  (*this) = Var<T>(vars,name);
  if (!exists()) (*this) = init;
}


template <class T>
Var<T>::Var(Variables *vars, const std::string &inname) :
  scope_(vars), scope_index_(0)
{
  std::string name = inname;
  if (vars->alias_.find(name) != vars->alias_.end()) {
    name = vars->alias_[name];
  }
  Variables::var_value_t varval = vars->find_value_ptr(name);

  Variables::var_type_e var_type = vars->type_to_enum(T());

  if (varval.second) {
    Variables::value_t *value_ptr = varval.second;
    if (value_ptr->var_type_ == Variables::UNKNOWN_E ||
        value_ptr->var_type_ == Variables::STRING_E ) {
      value_ptr->var_type_ = var_type;
    } else if (value_ptr->var_type_ != var_type) {
      throw "invalid type change";
    }

    if (value_ptr->cache_index_ == -1) {
      value_ptr->update_cache_from_string(varval.first);
    }
  } else {
    varval = vars->insert_variable(name,var_type,false);
  }
  scope_ = varval.first;
  scope_index_ = &varval.second->cache_index_;
}


} // end namespace Skinner
} // end namespace SCIRun

#endif // #define Skinner_Variables_H
