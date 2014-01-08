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

itcl::class SCIRun_Math_CreateScalarMatrix {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name CreateScalarMatrix
        set_defaults
    }

    method set_defaults {} {
      global $this-value 
      set $this-value  0.0
      global $this-init
      set $this-init 0
    }

    method maybeRestart { args } {
      if {[set $this-init] == 1} {
        $this-c needexecute
      }
    }


    method update_range { args } {
    
      set w .ui[modname]
      if {[winfo exists $w]} {
      
        set val [$w.value childsite]
        upvar \#0 $this-value-min min $this-value-max max $this-value-step step 
        
        $val.slider configure -from $min -to $max -resolution $step -variable $this-value
      }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w

        sci_labeledframe $w.value -labeltext "Value:"
        set val [$w.value childsite]

        set value [set $this-value]
        sci_scale $val.slider -variable $this-value -length 200 -showvalue true -orient horizontal  -from [set $this-value-min] -to [set $this-value-max] 
        sci_spinint $val.count -range {0 86400000} -justify right -width 5 -step 1 -textvariable $this-value -repeatdelay 300 -repeatinterval 10
        set $this-value $value
        
        pack $val.slider -side left -fill x -expand yes -anchor w 
        pack $val.count -side left -padx 10

        pack $w.value -side top -fill x -anchor w
        
        sci_labeledframe $w.settings -labeltext "Settings:"
        set settings [$w.settings childsite]
        
        sci_label $settings.lmin -text "Min Value:"
        sci_label $settings.lmax -text "Max Value:"
        sci_label $settings.lstep -text "Step:"
        
        sci_entryfield $settings.min -textvariable $this-value-min -command "$this update_range"
        sci_entryfield $settings.max -textvariable $this-value-max -command "$this update_range"
        sci_entryfield $settings.step -textvariable $this-value-step -command "$this update_range"
        
        $this update_range
                
        grid $settings.lmin -row 0 -column 0 -stick w
        grid $settings.lmax -row 1 -column 0 -stick w
        grid $settings.lstep -row 2 -column 0 -stick w

        grid $settings.min -row 0 -column 1 -stick w
        grid $settings.max -row 1 -column 1 -stick w
        grid $settings.step -row 2 -column 1 -stick w

        pack $w.settings -side top -anchor w -fill x
        
        $val.count configure -command "$this maybeRestart"
        $val.slider configure -command "$this maybeRestart"

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}


