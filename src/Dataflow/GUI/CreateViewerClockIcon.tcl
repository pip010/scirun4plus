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

#
#  CreateViewerClockIcon.tcl
#
#  Written by:
#   Allen Sanderson
#   SCI Institute
#   University of Utah
#   January 2004

itcl::class SCIRun_Visualization_CreateViewerClockIcon {
    inherit Module
   
     constructor { {args ""} } {
        eval configure $args
        set name CreateViewerClockIcon
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w
	wm title $w "Gen Clock"
	
# Type
	sci_labeledframe $w.type -labeltext "Clock Type"
	set type [$w.type childsite]

# Type - analog
	sci_frame $type.analog

	sci_radiobutton $type.analog.button -variable $this-type -value 0 \
	     -command "$this-c needexecute"
	sci_label $type.analog.label -text "Analog" -width 8 \
	    -anchor w -just left
	
	pack $type.analog.button $type.analog.label -side left

# Type - digital
	sci_frame $type.digital

	sci_radiobutton $type.digital.button -variable $this-type -value 1 \
	     -command "$this-c needexecute"
	sci_label $type.digital.label -text "Digital" -width 8 \
	    -anchor w -just left
	
	pack $type.digital.button $type.digital.label -side left

# Type - analog/digital
	sci_frame $type.both

	sci_radiobutton $type.both.button -variable $this-type -value 2 \
	     -command "$this-c needexecute"
	sci_label $type.both.label -text "Both" -width 8 \
	    -anchor w -just left
	
	pack $type.both.button $type.both.label -side left


	pack $type.analog $type.digital $type.both -side left
	pack $w.type -fill x -expand yes -side top


# Style
	sci_labeledframe $w.style -labeltext "Clock Style"
	set style [$w.style childsite]

# Style - box
	sci_frame $style.bbox

	sci_checkbutton $style.bbox.button -variable $this-bbox \
	     -command "$this-c needexecute"
	sci_label $style.bbox.label -text "Box" -width 4 \
	    -anchor w -just left
	
	pack $style.bbox.button $style.bbox.label -side left

	pack $style.bbox -side left

# Style - color
	sci_frame $style.color
	addColorSelection $style.color "Color" $this-color "color_change"
	pack $style.color -side left -padx 5

	pack $w.style -fill x -expand yes -side top

# Style - format
	sci_frame $style.format
	sci_label $style.format.label -text "C Style Format" -width 15 \
	    -anchor w -just left
	sci_entry $style.format.entry -width 16 -text $this-format

	pack $style.format.label $style.format.entry -side left

	pack $style.format -side left -padx 5


# Range
	sci_labeledframe $w.range -labeltext "Analog Clock Range"
	set range [$w.range childsite]

# Range - minimum
	sci_frame $range.min
	sci_label $range.min.label -text "Min."  -width 5 -anchor w -just left
	sci_entry $range.min.entry -width 6 -text $this-min

	pack $range.min.label $range.min.entry -side left
	pack $range.min -side left

# Range - maximum
	sci_frame $range.max
	sci_label $range.max.label -text "Max."  -width 5 -anchor w -just left
	sci_entry $range.max.entry -width 6 -text $this-max

	pack $range.max.label $range.max.entry -side left
	pack $range.max -side left -padx 5

# Range - current
	sci_frame $range.current
	sci_label $range.current.label -text "Current" \
	    -width 7 -anchor w -just left
	sci_entry $range.current.entry -width 6 -text $this-current

	pack $range.current.label $range.current.entry -side left
	pack $range.current -side left -padx 5

	pack $w.range -fill x -expand yes -side top


	sci_frame $w.twocol
# -relief groove -borderwidth 2
	
# Size
	sci_labeledframe $w.twocol.size -labeltext "Clock Size"
	set size [$w.twocol.size childsite]

# Size - tiny
        sci_frame $size.tiny

        sci_radiobutton $size.tiny.button -variable $this-size -value 0 \
            -command "$this-c needexecute"
        sci_label $size.tiny.label -text "Tiny" -width 6 \
            -anchor w -just left
        
        pack $size.tiny.button $size.tiny.label -side left
        pack $size.tiny -side top -padx 5

# Size - small
	sci_frame $size.small

	sci_radiobutton $size.small.button -variable $this-size -value 1 \
	    -command "$this-c needexecute"
	sci_label $size.small.label -text "Small" -width 6 \
	    -anchor w -just left
	
	pack $size.small.button $size.small.label -side left
	pack $size.small -side top -padx 5

# Size - medium
	sci_frame $size.medium

	sci_radiobutton $size.medium.button -variable $this-size -value 2 \
	    -command "$this-c needexecute"
	sci_label $size.medium.label -text "Medium" -width 6 \
	    -anchor w -just left
	
	pack $size.medium.button $size.medium.label -side left
	pack $size.medium -side top -padx 5

# Size - large
	sci_frame $size.large

	sci_radiobutton $size.large.button -variable $this-size -value 3 \
	    -command "$this-c needexecute"
	sci_label $size.large.label -text "Large" -width 6 \
	    -anchor w -just left
	
	pack $size.large.button $size.large.label -side left
	pack $size.large -side top -padx 5

# Size - huge
        sci_frame $size.huge

        sci_radiobutton $size.huge.button -variable $this-size -value 4 \
            -command "$this-c needexecute"
        sci_label $size.huge.label -text "Huge" -width 6 \
            -anchor w -just left
        
        pack $size.huge.button $size.huge.label -side left
        pack $size.huge -side top -padx 5

# Location
	sci_labeledframe $w.twocol.location -labeltext "Clock Location"
	set location [$w.twocol.location childsite]

	set locator [makeStickyLocator $location.gui \
			 $this-location-x $this-location-y \
			 100 100]

	$locator bind movable <ButtonRelease> "$this-c needexecute"

	pack $location.gui -fill x -expand yes -side top


	pack $w.twocol.size $w.twocol.location -fill both -expand yes -side left

	pack $w.twocol -fill x -expand yes -side top

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
	 frame $frame.colorFrame.col -relief ridge -borderwidth \
		 4 -height 0.8c -width 1.0c \
		 -background [format #%04x%04x%04x $ir $ig $ib]
	 
	 set cmmd "$this raiseColor $frame.colorFrame.col $color $colMsg"
	 sci_button $frame.colorFrame.set_color \
		 -text $text -command $cmmd
	 
	 #pack the node color frame
	 pack $frame.colorFrame.set_color $frame.colorFrame.col \
				  -side left -padx 2
	 pack $frame.colorFrame -side left
    }
}

