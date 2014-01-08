#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2012 Scientific Computing and Imaging Institute,
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


itcl::class SCIRun_Math_ReportRootMeanSquareError {
  inherit Module

  constructor { {args ""} } {
    eval configure $args
      set name ReportRootMeanSquareError
    }

    method matrix_order_widgets_enabled {} {
      set w .ui[modname]
      if {![winfo exists $w]} {
        return
      }

      $w.f.matrix_order_frame.row_order configure -state normal
      $w.f.matrix_order_frame.column_order configure -state normal
    }

    method matrix_order_widgets_disabled {} {
      set w .ui[modname]
      if {![winfo exists $w]} {
        return
      }

      $w.f.matrix_order_frame.row_order configure -state disabled
      $w.f.matrix_order_frame.column_order configure -state disabled
    }

    method ui {} {
      set w .ui[modname]
      if {[winfo exists $w]} {
        raise $w
        return
      }

      sci_toplevel $w
      wm minsize $w 160 30

      sci_frame $w.f
      pack $w.f -padx 2 -pady 2 -side top -expand yes

      sci_frame $w.f.matrix_order_frame
      pack $w.f.matrix_order_frame -padx 2 -pady 2 -side top -expand yes

      sci_frame $w.f.rmse_frame
      pack $w.f.rmse_frame -padx 2 -pady 2 -side top -expand yes

      sci_label $w.f.matrix_order_frame.label -text "Matrix Order: "
      pack $w.f.matrix_order_frame.label -padx 2 -pady 2 -side left

      # row order = 0, column order = 1 (C++ class)
      sci_radiobutton $w.f.matrix_order_frame.row_order -text "Row order" \
        -variable "$this-matrix_order" -value 0 -state disabled
      pack $w.f.matrix_order_frame.row_order -padx 2 -pady 2 -side right

      sci_radiobutton $w.f.matrix_order_frame.column_order -text "Column order" \
        -variable "$this-matrix_order" -value 1 -state disabled
      pack $w.f.matrix_order_frame.column_order -padx 2 -pady 2 -side right

     if {[set $this-input_is_vector] == 0} {
       matrix_order_widgets_enabled
     }

      sci_label $w.f.rmse_frame.rmse_label -text "Root Mean Square Error: "
      pack $w.f.rmse_frame.rmse_label -padx 2 -pady 2 -side left -expand yes

      sci_entry $w.f.rmse_frame.rmse_result -textvariable $this-rmse -state readonly
      pack $w.f.rmse_frame.rmse_result -padx 2 -pady 2 -side right -expand yes

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }
}
