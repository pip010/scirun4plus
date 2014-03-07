##
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

itcl::class SCIRun_Visualization_ShowString {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name ShowString
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w
 
        # Style
        sci_labeledframe $w.style -labeltext "String Style"
        set style [$w.style childsite]

        # Style - color
        sci_frame $style.color
        addColorSelection $style.color "Color" $this-color "color_change"
        pack $style.color -side left -padx 5

        pack $w.style  -fill x -expand yes -side top 

        sci_frame $w.twocol
        
        # Size
        sci_labeledframe $w.twocol.size -labeltext "String Size"
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
        sci_labeledframe $w.twocol.location -labeltext "String Location"
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
       pack $frame.colorFrame.set_color $frame.colorFrame.col -side left -padx 2
       pack $frame.colorFrame -side left
   }
}

