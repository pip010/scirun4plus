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

// http://www.acm.org/pubs/tog/GraphicsGems/gemsiv/arcball
// adapted from  Ken Shoemakes graphics gems code into c++, PPS
/* Ken Shoemake, 1993 */

#if !defined(Tools_Ball_h)
#define Tools_Ball_h

#include <Core/Events/OGLView/BallAux.h>

namespace SCIRun {

typedef enum AxisSet{NoAxes, CameraAxes, BodyAxes, OtherAxes, NSets} AxisSet;
typedef double *ConstraintSet;
class BallData {
public:
  BallData() { Init(); }
  void Init(void);
  void Place(HVect cntr, double rad) {
    center = cntr;
    radius = rad;
  }
  void Mouse(HVect vnow) { vNow = vnow; }
  void UseSet(AxisSet);
  void ShowResult(void);
  void HideResult(void);
  void Update(void);
  void Value(HMatrix&);
  void BeginDrag(void);
  void EndDrag(void);
  void Draw(void);
  void DValue(HMatrix&); // value at qDown
  void SetAngle(double); // works of normalized value...

  HVect center;
  double radius;
  BallQuaternion qNow, qDown, qDrag, qNorm;
  HVect vNow, vDown, vFrom, vTo, vrFrom, vrTo;
  HMatrix mNow, mDown;
  int showResult, dragging;
  ConstraintSet sets[NSets];
  int setSizes[NSets];
  AxisSet axisSet;
  int axisIndex;

private:
  static void DrawAnyArc(HVect vFrom, HVect vTo);
  static void DrawHalfArc(HVect n);
  static void Ball_DrawConstraints(BallData *ball);
  static void Ball_DrawDragArc(BallData *ball);
  static void Ball_DrawResultArc(BallData *ball);
};

} // End namespace SCIRun




#endif

