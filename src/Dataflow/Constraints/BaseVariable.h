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
 *  BaseVariable.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Aug. 1994
 *
 */


#ifndef SCI_project_Base_Variable_h
#define SCI_project_Base_Variable_h 1

#include <Dataflow/Constraints/VarCore.h>
#include <string>
#include <vector>

#include <Dataflow/Constraints/share.h>

namespace SCIRun {

/* Priority levels */
// P_constant is for a variable in reference to a constraint.
enum VPriority { P_Lowest, P_LowMedium, P_Default, P_HighMedium, P_Highest };

enum Scheme { Scheme1, Scheme2, Scheme3, Scheme4,
	      Scheme5, Scheme6, Scheme7, Scheme8, DefaultScheme };

// USE THESE TO BE CLEAR!
#define PointVariable BaseVariable
#define RealVariable BaseVariable

class BaseConstraint;
class ConstraintSolver;

class SCISHARE BaseVariable {
   friend class BaseConstraint;
   friend class ConstraintSolver;
public:
   BaseVariable( const std::string& name, ConstraintSolver* s, const Scheme scheme,
		 const Point& initialValue );
   BaseVariable( const std::string& name, ConstraintSolver* s, const Scheme scheme,
		 const double initialValue );
   ~BaseVariable();

   void Order(); // Use to let the Variable order its constraints.

   // Widgets use these instead of Assign!
   void Set( const Point& newValue, const Scheme s = DefaultScheme );
   void Set( const double newValue, const Scheme s = DefaultScheme );
   void SetDelta( const Vector& deltaValue, const Scheme s = DefaultScheme );
   void SetDelta( const double  deltaValue, const Scheme s = DefaultScheme );

   // Widgets use these to move whole widget.
   // i.e. they don't change constraints!!
   void Move( const Point& newValue );
   void Move( const double newValue );
   void MoveDelta( const Vector& deltaValue );
   void MoveDelta( const double deltaValue );
   
   inline Point point() const;
   inline operator Point() const;
   inline double real() const;
   inline operator double() const;

   void print( std::ostream& os );

   inline const std::string& GetName() const { return name; }
   inline int GetNumConstraints() const { return numconstraints; }

private:
   std::string name;

   VarCore data;

   ConstraintSolver* solver;
   
   unsigned int level;

   Index numconstraints;
   std::vector<BaseConstraint*> constraints;
   std::vector<Index>           constraint_indices;
   std::vector<VPriority>       constraint_priorities;
   std::vector<Index>           constraint_order;
   
   Scheme scheme;

   Index Register( BaseConstraint* constraint, const Index index );
   void RegisterPriority( const Index index, const VPriority p );
   void printc( std::ostream& os, const Index c );
};
inline std::ostream& operator<<( std::ostream& os, BaseVariable& v );


/* Miscellaneous */
const char* PriorityString( const VPriority p );
const char* SchemeString( const Scheme s );

inline int HigherPriority( const VPriority p1, const VPriority p2 )
{
   return (p1 > p2);
}


inline Point
BaseVariable::point() const
{
  return static_cast<Point>(data);
}


inline
BaseVariable::operator Point() const
{
  return static_cast<Point>(data);
}


inline double
BaseVariable::real() const
{
  return static_cast<double>(data);
}


inline
BaseVariable::operator double() const
{
  return static_cast<double>(data);
}


inline std::ostream&
operator<<( std::ostream& os, BaseVariable& v )
{
  v.print(os);
  return os;
}

} // End namespace SCIRun


#endif
