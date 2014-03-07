##
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

itcl::class SCIRun_ChangeFieldData_RegisterWithCorrespondences {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name RegisterWithCorrespondences
        set_defaults
    }

    method set_defaults {} {
    	global $this-method
        global $this-ed-method
        
        set $this-method "transform"
        set $this-ed-method "affine"
        
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w
        
        sci_frame $w.frame
        pack $w.frame -fill x
        
        sci_tabnotebook $w.frame.filters -height 100 -tabpos n        
        $w.frame.filters add -label "Transform" -command "set $this-method transform"
        $w.frame.filters select 0

        pack $w.frame.filters -fill x -expand yes
        
        set transform  [$w.frame.filters childsite 0]
 

        sci_frame $transform.f
        pack $transform.f -side left -padx 5p -anchor n

        sci_radiobutton $transform.f.affine -text "Affine" -variable $this-ed-method -value "affine"
        sci_radiobutton $transform.f.morph -text "Morph" -variable $this-ed-method -value "morph"
        sci_radiobutton $transform.f.none -text "None" -variable $this-ed-method -value "none"
     
        
        grid $transform.f.affine -row 0 -column 0 -sticky w
        grid $transform.f.morph -row 0 -column 1 -sticky w
        grid $transform.f.none -row 1 -column 0 -columnspan 2 -sticky w
      
        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}


