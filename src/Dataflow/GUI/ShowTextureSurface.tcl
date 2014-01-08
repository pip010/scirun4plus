#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#


catch {rename ShowTextureSurface ""}

itcl::class SCIRun_Visualization_ShowTextureSurface {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name ShowTextureSurface
    }

    method ui {} { 
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }
        sci_toplevel $w

        set n "$this-c needexecute "

        sci_frame $w.main -relief flat
        pack $w.main -fill both -expand yes

        sci_frame $w.main.f3 -relief groove -borderwidth 2
        pack $w.main.f3 -padx 4 -pady 4 -fill x -side top

        sci_label $w.main.f3.l -text "Interpolation Mode"
        sci_radiobutton $w.main.f3.interp -text "Interpolate" -relief flat \
          -variable $this-interp_mode -value 1 \
          -anchor w -command $n

        sci_radiobutton $w.main.f3.near -text "Nearest" -relief flat \
          -variable $this-interp_mode -value 0 \
          -anchor w -command $n

        pack $w.main.f3.l $w.main.f3.interp $w.main.f3.near -side top -fill x -padx 4 -pady 2
        
        makeSciButtonPanel $w.main $w $this
        moveToCursor $w
    }
}

