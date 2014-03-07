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

catch {rename SCIRun_NewField_RefineMesh ""}

itcl::class SCIRun_NewField_RefineMesh {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name RefineMesh
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return;
        }

        sci_toplevel $w
        #wm minsize $w 80 130

        sci_frame $w.style -relief groove -borderwidth 2
        sci_label $w.style.t1 -text "Refinement Method"
        pack $w.style.t1

        sci_radiobutton $w.style.default -text "Default" \
            -variable $this-method -value "default"
        sci_radiobutton $w.style.convex -text "Expand refinement volume to improve element quality" \
            -variable $this-method -value "convex"
        pack $w.style.default $w.style.convex \
            -side top -anchor w

        sci_frame $w.isovalue -relief groove -borderwidth 2
        sci_label $w.isovalue.t1 -text "Addional Constraints"
        pack $w.isovalue.t1 

        sci_radiobutton $w.isovalue.all -text "Do not add constraint" \
            -variable $this-select -value "all"
        sci_radiobutton $w.isovalue.lessthan -text "Do not refine nodes/elements with values less than isovalue" \
            -variable $this-select -value "greaterthan"
        sci_radiobutton $w.isovalue.equalto -text "Do not refine nodes/elements with values unequal to isovalue" \
            -variable $this-select -value "equal"
        sci_radiobutton $w.isovalue.greaterthan -text "Do not refine nodes/elements with values greater than isovalue" \
            -variable $this-select -value "lessthan"
        sci_radiobutton $w.isovalue.none -text "Do not refine any elements" \
            -variable $this-select -value "none"

        sci_frame $w.isovalue.f 
        
        sci_label $w.isovalue.f.label -text "IsoValue :"
        sci_entry $w.isovalue.f.value  -textvariable $this-isoval
        pack $w.isovalue.f.label $w.isovalue.f.value -side left -anchor n
                
        pack $w.isovalue.all $w.isovalue.lessthan $w.isovalue.equalto $w.isovalue.greaterthan $w.isovalue.none $w.isovalue.f\
            -side top -anchor w

        pack $w.style $w.isovalue -side top -e y -f both -padx 5 -pady 5

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
