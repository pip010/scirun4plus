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
 *  FontManager.h
 *
 *  Written by:
 *   McKay Davis
 *   School of Computing
 *   University of Utah
 *   April, 2006
 *
 */

#ifndef SCIRun_Dataflow_Modules_Render_FontManager_h
#define SCIRun_Dataflow_Modules_Render_FontManager_h

#include <Core/Geom/FreeType.h>
#include <Core/Geom/TextRenderer.h>
#include <string>
#include <map>
#include <set>

#include <Core/Geom/share.h>
namespace SCIRun {

class SCISHARE FontManager {
  public:
    static TextRenderer * get_renderer(double, const std::string& filename = "scirun.ttf");
    //! For backwards compatibility, take a fontindex
    //! and choose between 5 fixed font sizes
    static TextRenderer * get_sized_renderer(int, const std::string& filename="scirun.ttf");

    static void   release_renderer(double, const std::string& filename="scirun.ttf");
    static void   release_renderer(TextRenderer *);
    static void		release_all();

  private:
    
    static FreeTypeFace*         load_face(const std::string &filename);
    // Pointer to the freetype library
    static FreeTypeLibraryHandle   freetype_lib_;

    typedef  std::map<int, TextRenderer *> SizeRendererMap_t;
    typedef  std::map<std::string, SizeRendererMap_t> NameSizeRendererMap_t;
    typedef  std::map<TextRenderer *, FreeTypeFaceHandle> RendererFaceMap_t;
    
    static   NameSizeRendererMap_t  renderers_;
    static   RendererFaceMap_t      face_;
};

}


#endif
