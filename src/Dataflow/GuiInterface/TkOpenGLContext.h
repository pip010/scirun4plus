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
 *  TkOpenGLContext.h:
 *
 *  Written by:
 *   McKay Davis
 *   December 2004
 */


#ifndef SCIRun_Core_2d_TkOpenGLContext_h
#define SCIRun_Core_2d_TkOpenGLContext_h

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Core/Geom/OpenGLContext.h>
#include <Dataflow/GuiInterface/share.h>

namespace SCIRun {

class TkOpenGLContextPrivate;

class SCISHARE TkOpenGLContext : public OpenGLContext 
{
  public:
    TkOpenGLContext(const std::string &, int visualid=0, 
        int width=640, int height = 480);
    virtual ~TkOpenGLContext();

    std::string       listvisuals();
    virtual bool			make_current();
    virtual void			release();
    virtual int       width();
    virtual int       height();
    virtual void			swap();
    virtual void			lock();
    virtual void			unlock();
    virtual bool      has_shaders();
    virtual bool      initialized();
    void restackWindow();
    void define_cursor(const std::string& name);
    void xsync();
    boost::shared_ptr<TkOpenGLContextPrivate> impl_;
  private:
    friend class TkOpenGLContextPrivate;
    std::vector<int>	valid_visuals_;
    int               visualid_;
    int               screen_number_;
    std::string       id_;
};

} // End namespace SCIRun

#endif // SCIRun_Core_2d_OpenGLContext_h
