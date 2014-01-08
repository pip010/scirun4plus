#
#  For more information, please see: http://software.sci.utah.edu
# 
# The MIT License
#
# Copyright (c) 2009 Scientific Computing and Imaging Institute,
# University of Utah.
#
#
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
#

itcl::class SCIRun_FieldArray_ExtractFieldArrayRange {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name ExtractFieldArrayRange
        set_defaults
    }

    method set_defaults {} {
    }

    method update_range { } {
    
      set w .ui[modname]
      if {[winfo exists $w]} {
      
        upvar \#0 $this-range-min min $this-range-max max         
        $w.f1.start configure -from $min -to $max -variable $this-range-start
        $w.f1.starts configure -command "$this update_range "
        $w.f2.end configure -from $min -to $max -variable $this-range-end 
        $w.f2.ends configure -command "$this update_range "
      }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        toplevel $w

        sci_frame $w.f1
        sci_frame $w.f2
        pack $w.f1 $w.f2 -side top -anchor n -fill x -expand yes
        
        set tmp [set $this-range-start]
        sci_scale $w.f1.start -from 0 -to [set $this-range-end]  -variable $this-range-start -showvalue true -orient horizontal -relief groove -label "Field Index Start:" 
        sci_spinint $w.f1.starts -range {0 86400000} -justify right -width 5 -step 1 -textvariable $this-range-start -repeatdelay 300 -repeatinterval 10 
        pack $w.f1.start -side left -anchor n -fill x -expand yes -padx 8 -pady 5
        pack $w.f1.starts -side left -padx 8 -pady 5
        set $this-range-start $tmp

        set tmp [set $this-range-end]
        sci_scale $w.f2.end -from 0 -to [set $this-range-end] -variable $this-range-end -showvalue true -orient horizontal -relief groove -label "Field Index End:" 
        sci_spinint $w.f2.ends -range {0 86400000} -justify right -width 5 -step 1 -textvariable $this-range-end -repeatdelay 300 -repeatinterval 10 
        pack $w.f2.end -side left -anchor n -fill x -expand yes -padx 8 -pady 5
        pack $w.f2.ends -side left -padx 8 -pady 5
        set $this-range-end $tmp

        update_range

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}


