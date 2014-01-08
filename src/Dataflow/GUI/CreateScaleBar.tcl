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


itcl::class SCIRun_Visualization_CreateScaleBar {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name CreateScaleBar
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        # Text
        sci_frame $w.text -relief groove -borderwidth 2

        # Style - text color
        frame $w.text.color
        addColorSelection $w.text.color "Text Color" $this-color "color_change"

        pack $w.text.color -side top -fill x -pady 3

        # Style - Labeling
        sci_frame $w.tc_uns -relief groove -borderwidth 2
        sci_frame $w.tc_uns.left
        sci_frame $w.tc_uns.right

        sci_frame $w.tc_uns.left.units -borderwidth 2
        sci_label $w.tc_uns.left.units.label -text "Label"
        sci_entry $w.tc_uns.left.units.entry -width 10 -textvariable $this-label
        bind $w.tc_uns.left.units.entry <KeyPress-Return> "$this-c needexecute"

        sci_frame $w.tc_uns.right.scale -borderwidth 2
        sci_label $w.tc_uns.right.scale.label -text "Scale"
        sci_entry $w.tc_uns.right.scale.entry -width 10 -textvariable $this-scale
        bind $w.tc_uns.right.scale.entry <KeyPress-Return> "$this-c needexecute"

        # Pack the Labels/Units/Scale Widgets Frame
        pack $w.tc_uns.left.units.label   $w.tc_uns.left.units.entry   -side left -anchor e
        pack $w.tc_uns.right.scale.label   $w.tc_uns.right.scale.entry   -side left -anchor e

        pack $w.tc_uns.left.units \
            -side top -padx 5 -pady 2 -anchor e

        pack $w.tc_uns.right.scale \
            -side top -padx 5 -pady 2 -anchor e

        pack $w.tc_uns.left $w.tc_uns.right -side left

        pack $w.text $w.tc_uns \
            -anchor w -side top -fill x -pady 2

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }

    method raiseColor {col color colMsg} {
        global $color
        set window .ui[modname]
        if {[winfo exists $window.color]} {
            SciRaise $window.color
            return
        } else {
            makeColorPicker $window.color $color \
          "$this setColor $col $color $colMsg" \
          "destroy $window.color"
        }
    }

    method setColor {col color colMsg} {
        global $color
        global $color-r
        global $color-g
        global $color-b
        set ir [expr int([set $color-r] * 65535)]
        set ig [expr int([set $color-g] * 65535)]
        set ib [expr int([set $color-b] * 65535)]
        
        set window .ui[modname]
        $col config -background [format #%04x%04x%04x $ir $ig $ib]
        $this-c $colMsg
               
              # The above works for only the geometry not for the text so execute.
        $this-c needexecute
    }

    method addColorSelection {frame text color colMsg} {
        #add node color picking 
        global $color
        global $color-r
        global $color-g
        global $color-b
        set ir [expr int([set $color-r] * 65535)]
        set ig [expr int([set $color-g] * 65535)]
        set ib [expr int([set $color-b] * 65535)]
        
        sci_frame $frame.colorFrame
        sci_frame $frame.colorFrame.col -relief ridge -borderwidth \
            4 -height 0.8c -width 1.0c \
            -background [format #%04x%04x%04x $ir $ig $ib]
             
        set cmmd "$this raiseColor $frame.colorFrame.col $color $colMsg"
             button $frame.colorFrame.set_color \
                 -text $text -command $cmmd
             
        #pack the node color frame
        pack $frame.colorFrame.set_color $frame.colorFrame.col -side left -padx 2
        pack $frame.colorFrame -side left
    }
}
  



