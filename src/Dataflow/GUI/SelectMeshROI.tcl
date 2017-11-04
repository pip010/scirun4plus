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

catch {rename SCIRun_MiscField_SelectMeshROI ""}

itcl::class SCIRun_MiscField_SelectMeshROI {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name SelectMeshROI
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return;
        }

        sci_toplevel $w
        #wm minsize $w 80 130


        sci_frame $w.isovalue -relief groove -borderwidth 2
        sci_label $w.isovalue.t1 -text "Select ROI Method"
        pack $w.isovalue.t1 

        sci_radiobutton $w.isovalue.point -text "Point in space" \
            -variable $this-select -value "point"
        sci_radiobutton $w.isovalue.line -text "Line intersection" \
            -variable $this-select -value "line"

        sci_frame $w.isovalue.f 
        
        sci_label $w.isovalue.f.label -text "Topological distance (range) :"
        sci_entry $w.isovalue.f.value -textvariable $this-isoval
        pack $w.isovalue.f.label $w.isovalue.f.value -side left -anchor n
                
        pack $w.isovalue.point $w.isovalue.line $w.isovalue.f -side top -anchor w

        pack $w.isovalue -side top -e y -f both -padx 5 -pady 5

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
