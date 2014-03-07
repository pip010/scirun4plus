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
 *  GuiGeom.h: Interface to TCL variables for Geom stuff
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   June 1995
 *
 */

#ifndef SCI_Geom_GuiGeom_h
#define SCI_Geom_GuiGeom_h 1

#include <Core/Datatypes/Color.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/GuiInterface/share.h>

namespace SCIRun {

class SCISHARE GuiColor : public GuiVar {
  private:
    GuiDouble r_;
    GuiDouble g_;
    GuiDouble b_;
  public:
    GuiColor(GuiContext* ctx);
    GuiColor(GuiContext* ctx,Color c);
    
    ~GuiColor();

    Color get();
    Color get_cached();
    void  request();
    
    void set(const Color&);
    void send(const Color&);
  };

class Material;
  
class SCISHARE GuiMaterial : public GuiVar {
  private:
    GuiColor  ambient_;
    GuiColor  diffuse_;
    GuiColor  specular_;
    GuiDouble shininess_;
    GuiColor  emission_;
    GuiDouble reflectivity_;
    GuiDouble transparency_;
    GuiDouble refraction_index_;
    
  public:
    GuiMaterial(GuiContext* ctx);
    ~GuiMaterial();
   
    Material get();
    Material get_cached();
    void     request();
    
    void set(const Material&);
    void send(const Material&);
};

} // End namespace SCIRun


#endif // ifndef SCI_Geom_GuiGeom_h
