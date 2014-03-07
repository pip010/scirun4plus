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

itcl::class SCIRun_ChangeMesh_PadRegularMesh {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name PadRegularMesh
        set_defaults
    }

    method set_defaults {} {
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w
        
        sci_frame $w.pad 
        pack $w.pad -fill x -expand yes -side top

        sci_label $w.pad.lsize -text "Number of elements to be padded:"
        sci_label $w.pad.lvalue -text "Value inserted into padded region:"

        sci_entry $w.pad.size -textvariable $this-pad-size
        sci_entry $w.pad.value -textvariable $this-pad-value
        
        grid $w.pad.lsize -column 0 -row 0 -sticky w
        grid $w.pad.lvalue -column 0 -row 1 -sticky w
        grid $w.pad.size -column 1 -row 0 -sticky w
        grid $w.pad.value -column 1 -row 1 -sticky w

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}


