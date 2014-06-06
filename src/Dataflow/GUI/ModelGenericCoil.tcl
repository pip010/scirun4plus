#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2014 Scientific Computing and Imaging Institute,
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

itcl::class SCIRun_Math_ModelGenericCoil {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
      set name ModelGenericCoil
      set_defaults
    }

    method set_defaults {} {
      global $this-wireCurrentTCL
      global $this-coilRadiusTCL
      global $this-coilDistanceTCL
      global $this-coilSegmentsTCL
      global $this-typeTCL
      set $this-wireCurrentTCL 1
      set $this-coilRadiusTCL 10
	  set $this-coilDistanceTCL 2
	  set $this-coilSegmentsTCL 33
      set $this-typeTCL "O-shaped"
    }
 
    method make_entry {w text v c} {
        sci_frame $w
        sci_label $w.l -text "$text"
        pack $w.l -side left
        sci_entry $w.e -textvariable $v
        bind $w.e <Return> $c
        pack $w.e -side right
    }
    method ui {} {
      global $this-wireCurrentTCL
      global $this-coilRadiusTCL
	  global $this-coilDistanceTCL
      global $this-coilSegmentsTCL
      global $this-typeTCL

      set w .ui[modname]
      if {[winfo exists $w]} {
          return
      }
      sci_toplevel $w
      
      make_labeled_radio $w.type "Coil type:" "" left 2 $this-typeTCL \
          {"O-shaped" {"8-shaped"}}
      make_entry $w.current "Current through wire:" $this-wireCurrentTCL "$this-c needexecute"
      make_entry $w.radius "Radius of coil(s):" $this-coilRadiusTCL "$this-c needexecute"
      make_entry $w.distance "Distance between centers:" $this-coilDistanceTCL "$this-c needexecute"
      make_entry $w.segments "Coil segments:" $this-coilSegmentsTCL "$this-c needexecute"
      
      bind $w.current <Return> "$this-c needexecute"
      bind $w.radius <Return> "$this-c needexecute"
      bind $w.distance <Return> "$this-c needexecute"
      bind $w.segments <Return> "$this-c needexecute"
      
      pack $w.type $w.current $w.radius $w.distance $w.segments -side top -fill x

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }
 
    


}
