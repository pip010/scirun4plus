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

itcl::class SCIRun_TMS_ModelTMSCoilSingle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
      set name ModelTMSCoilSingle
      set_defaults
    }

    method set_defaults {} {
      global $this-wireCurrentTCL
      global $this-coilRadiusTCL
      global $this-levelDetailTCL
      global $this-outerDistanceTCL
      global $this-coilLayersTCL
      global $this-typeTCL
      
      set $this-wireCurrentTCL 1
      set $this-coilRadiusTCL 0.035
	    set $this-levelDetailTCL 6
      set $this-outerDistanceTCL 0.002
      set $this-coilLayersTCL 1
      set $this-typeTCL "8-shape"
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
      global $this-levelDetailTCL
      global $this-outerDistanceTCL
      global $this-coilLayersTCL
      global $this-typeTCL

      set w .ui[modname]
      if {[winfo exists $w]} {
          return
      }
      sci_toplevel $w
      
      make_labeled_radio $w.type "Type :" "" left 2 $this-typeTCL \
          { "0-shape" "8-shape" }
      make_entry $w.current "Current:" $this-wireCurrentTCL "$this-c needexecute"
      make_entry $w.radius "Radius :" $this-coilRadiusTCL "$this-c needexecute"
      make_entry $w.distance "Distance:" $this-outerDistanceTCL "$this-c needexecute"
      make_entry $w.layers "Layers:" $this-coilLayersTCL "$this-c needexecute"
      make_entry $w.lod "LOD:" $this-levelDetailTCL "$this-c needexecute"
      
      bind $w.current <Return> "$this-c needexecute"
      bind $w.radius <Return> "$this-c needexecute"
      bind $w.distance <Return> "$this-c needexecute"
      bind $w.layers <Return> "$this-c needexecute"
      bind $w.lod <Return> "$this-c needexecute"
      
      pack $w.type $w.current $w.radius $w.distance $w.layers $w.lod -side top -fill x

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }
 
    


}
