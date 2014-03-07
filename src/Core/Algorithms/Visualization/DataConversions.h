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
 *  DataConversions.cc: Data converters for rendering alogrithms
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   April 2007
 *
 */

#if !defined(Visualization_DataConversions_h)
#define Visualization_DataConversions_h

#include <Core/Datatypes/Color.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Tensor.h>

#include <sstream>
#include <iomanip>
#include <iostream>

#undef SCISHARE
#if defined(_WIN32) && !defined(BUILD_SCIRUN_STATIC)
#define SCISHARE __declspec(dllexport)
#else
#define SCISHARE
#endif

namespace SCIRun {

// Conversion templates.
template <class T>
SCISHARE bool
to_color( const T &v, Color &c )
{
  c = Color(fabs(v), fabs(v), fabs(v));
  return true;
}

template <class T>
SCISHARE bool
to_double( const T& data_in, double &data_out)
{
  data_out = static_cast<double>(data_in);
  return true;
}


template <class T>
SCISHARE bool
to_vector( const T&, Vector&)
{
  return false;
}


template <class T>
SCISHARE bool
to_tensor( const T&, Tensor&)
{
  return false;
}

template <class T>
SCISHARE bool
to_buffer( const T &value, std::ostringstream &buffer)
{
  buffer << value;
  return true;
}


// Conversion template specialization.
template <>
SCISHARE bool to_color( const Vector&, Color& );

template <>
SCISHARE bool to_color( const Tensor&, Color& );

template <>
SCISHARE bool to_double(const Vector&, double&);

template <>
SCISHARE bool to_double(const Tensor&, double&);

template <>
SCISHARE bool to_double(const std::string&, double&);

template <>
SCISHARE bool to_vector(const Vector&, Vector&);

template <>
SCISHARE bool to_vector(const Tensor&, Vector&);

template <>
SCISHARE bool to_tensor(const Tensor&, Tensor&);

template <>
SCISHARE bool to_buffer(const unsigned char&, std::ostringstream&);

template <>
SCISHARE bool to_buffer(const char&, std::ostringstream&);

}

#endif
