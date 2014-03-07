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

itcl::class SCIRun_NewField_GetDomainBoundary {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name GetDomainBoundary
        set_defaults
    }

    method set_defaults {} {
        global $this-userange
        global $this-usevalue
        global $this-minrange
        global $this-maxrange
        global $this-value
        global $this-includeouterboundary
        global $this-innerboundaryonly
        global $this-noinnerboundary
        global $this-disconnect
  
        set $this-userange 0
        set $this-usevalue 0
        set $this-minrange 0.0
        set $this-maxrange 255.0
        set $this-value 1.0
        set $this-includeouterboundary 1
        set $this-innerboundaryonly 0
        set $this-noinnerboundary 0
        set $this-disconnect 0
        
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w

        sci_frame $w.f
        sci_checkbutton $w.f.userange -text "Only include compartments in the range between:" -variable $this-userange -command "if { [set $this-usevalue] == 1 } { set $this-usevalue 0} ; $w.f.usevalue deselect"
        sci_checkbutton $w.f.usevalue -text "Only include compartment:" -variable $this-usevalue -command "if { [set $this-userange] == 1 } { set $this-userange 0} ; $w.f.userange deselect"

        sci_label $w.f.valuelabel -text "value:"
        sci_entry $w.f.value  -textvariable  $this-value

        sci_label $w.f.minrangelabel -text "min:"
        sci_entry $w.f.minrange  -textvariable  $this-minrange
        sci_label $w.f.maxrangelabel -text "max:"
        sci_entry $w.f.maxrange  -textvariable  $this-maxrange
        sci_checkbutton $w.f.includeouterboundary -text "Include outer boundary" -variable $this-includeouterboundary
        sci_checkbutton $w.f.innerboundaryonly -text "Include inner boundary only" -variable $this-innerboundaryonly
        sci_checkbutton $w.f.noinnerboundary -text "Exclude inner boundary" -variable $this-noinnerboundary
        sci_checkbutton $w.f.disconnect -text "Disconnect boundaries between different element types" -variable $this-disconnect

        grid $w.f.userange -column 0 -row 0 -columnspan 4 -sticky w
        grid $w.f.minrangelabel -column 0 -row 1 -sticky news
        grid $w.f.minrange -column 1 -row 1 -sticky news
        grid $w.f.maxrangelabel -column 2 -row 1 -sticky news
        grid $w.f.maxrange -column 3 -row 1 -sticky news
        grid $w.f.usevalue -column 0 -row 2 -columnspan 4 -sticky w
        grid $w.f.valuelabel -column 0 -row 3 -sticky news
        grid $w.f.value -column 1 -row 3 -sticky news
        grid $w.f.includeouterboundary -column 0 -row 4 -columnspan 4 -sticky w
        grid $w.f.innerboundaryonly -column 0 -row 5 -columnspan 4 -sticky w
        grid $w.f.noinnerboundary -column 0 -row 6 -columnspan 4 -sticky w
        grid $w.f.disconnect -column 0 -row 7 -columnspan 4 -sticky w

        pack $w.f -fill x
        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}


