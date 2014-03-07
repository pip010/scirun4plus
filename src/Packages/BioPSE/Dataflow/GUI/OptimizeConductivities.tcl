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


catch {rename BioPSE_Inverse_OptimizeConductivities ""}

itcl::class BioPSE_Inverse_OptimizeConductivities {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name OptimizeConductivities
        set_defaults
    }

    method set_defaults {} {	
        global $this-seed_gui
        set $this-seed_gui 0
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
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w
        wm minsize $w 100 50
        sci_frame $w.g
        sci_button $w.g.go -text "Execute" -relief raised -command "$this-c exec"
        sci_button $w.g.p -text "Pause" -relief raised -command "$this-c pause"
        sci_button $w.g.np -text "Unpause" -relief raised -command \
          "$this-c unpause"
        sci_button $w.g.stop -text "Stop" -relief raised -command "$this-c stop"
        pack $w.g.go $w.g.p $w.g.np $w.g.stop -side left -fill x -expand 1
        sci_frame $w.seed
        global $this-seed_gui
        make_entry $w.seed.seed "Random seed:" $this-seed_gui \
          "$this-c needexecute"
        pack $w.seed.seed -side top -fill x -expand 1
        pack $w.g $w.seed -side top -fill x -expand 1

        makeSciButtonPanel $w $w $this -no_execute
        moveToCursor $w
    }
}

