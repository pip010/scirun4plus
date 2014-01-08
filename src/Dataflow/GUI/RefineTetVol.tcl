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

#    File   : RefineTetVol.tcl
#    Author : Martin Cole
#    Date   : Thu Nov 13 10:10:00 2003

itcl::class SCIRun_FieldsCreate_RefineTetVol {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name RefineTetVol
    }

    method ui {} {
      set w .ui[modname]
      if {[winfo exists $w]} {
          raise $w
          return
      }
      sci_toplevel $w

      sci_frame $w.f

      make_labeled_radio $w.f.rb \
          "Execution mode" "" top 1 \
          $this-execution_mode \
          {{"Subdivide to level:" sub_all} {"Subdivide cell index:" sub_one}}


      sci_label $w.f.lab -text "Enter cell index or level"
      sci_entry $w.f.cell_index -width 40 -textvariable $this-cell_index

     	pack $w.f.rb $w.f.lab $w.f.cell_index -side top -anchor w
	
      sci_frame $w.controls
      sci_button $w.controls.exc -text "Execute" -command "$this-c needexecute"
      sci_button $w.controls.cancel -text "Cancel" -command "destroy $w"
      pack $w.controls.exc $w.controls.cancel -side left -fill both
      
      pack $w.f $w.controls -side top -expand yes -fill both -padx 5 -pady 5
    }
}




