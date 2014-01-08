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

##
 #  InsertHexVolSheetFromTriSurf.tcl: The InsertHexVolSheetAlongSurface UI
 #  Written by:
 #   Jason Shepherd
 #   Department of Computer Science
 #   University of Utah
 #   April 2006
 ##

itcl::class SCIRun_NewField_InsertHexVolSheetAlongSurface {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name InsertHexVolSheetAlongSurface
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return;
        }

        sci_toplevel $w
        wm minsize $w 150 80

        sci_frame $w.bound1
        sci_label $w.bound1.t1 -text "Intersected Hexes in"
        pack $w.bound1.t1
        pack $w.bound1

        sci_frame $w.bound
        sci_radiobutton $w.bound.side_1 -text "Side 1" \
            -variable $this-side -value "side1"
        sci_radiobutton $w.bound.side_2 -text "Side 2" \
            -variable $this-side -value "side2"
        pack $w.bound.side_1 $w.bound.side_2 \
            -side left -anchor nw -padx 3
        pack $w.bound -side top

        sci_frame $w.layer1
        sci_label $w.layer1.t1
        sci_label $w.layer1.t2 -text "Add Sheet?"
        pack $w.layer1.t1 $w.layer1.t2
        pack $w.layer1

        sci_frame $w.layer
        sci_radiobutton $w.layer.addlayeron -text "On" \
            -variable $this-addlayer -value "On"
        sci_radiobutton $w.layer.addlayeroff -text "Off" \
            -variable $this-addlayer -value "Off"
        pack $w.layer.addlayeron $w.layer.addlayeroff \
            -side left -anchor nw -padx 3
        pack $w.layer -side top

        sci_frame $w.f
        sci_frame $w.fb
        pack $w.f $w.fb -padx 2 -pady 2 -side top -expand yes

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
