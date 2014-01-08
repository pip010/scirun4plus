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


itcl::class SCIRun_MiscField_BuildMappingMatrix {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name BuildMappingMatrix
    }

    method ui {} {
      set w .ui[modname]
      if {[winfo exists $w]} {
          return
      }
      sci_toplevel $w

      sci_frame $w.basis
      sci_label $w.basis.label -text "Interpolation Basis:"
      sci_radiobutton $w.basis.const -text "Constant ('find closest')" \
        -variable $this-interpolation_basis -value constant
      sci_frame $w.basis.cframe 
      sci_label $w.basis.cframe.label -text "Constant Mapping:"
      sci_radiobutton $w.basis.cframe.onetomany -text \
        "Each destination gets nearest source value" \
        -variable $this-map_source_to_single_dest -value 0
      sci_radiobutton $w.basis.cframe.onetoone -text \
        "Each source projects to just one destination" \
        -variable $this-map_source_to_single_dest -value 1
      pack $w.basis.cframe.label -side top -anchor w
      pack $w.basis.cframe.onetomany $w.basis.cframe.onetoone \
        -side top -anchor w -padx 15
      sci_radiobutton $w.basis.lin -text "Linear (`weighted')" \
        -variable $this-interpolation_basis -value linear
      pack $w.basis.label -side top -anchor w
      pack $w.basis.const -padx 15 -side top -anchor w
      pack $w.basis.cframe -padx 30 -side top -anchor w
      pack $w.basis.lin -padx 15 -side top -anchor w
      
      sci_frame $w.exhaustive
      sci_label $w.exhaustive.label -text "Search Options:"
      sci_frame $w.exhaustive.dist
      sci_label $w.exhaustive.dist.label -text \
        "Maximum Distance (negative value -> 'no max'):"
      sci_entry $w.exhaustive.dist.entry -textvariable \
          $this-exhaustive_search_max_dist -width 8
      pack $w.exhaustive.dist.label $w.exhaustive.dist.entry \
          -side left -anchor n
      pack $w.exhaustive.label -side top -anchor w
      pack $w.exhaustive.dist -side top -anchor w -padx 30
      
      pack $w.basis -side top -anchor w
      pack $w.exhaustive -side top -anchor w -pady 15

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }
}
