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


itcl::class SCIRun_NewField_GenerateSinglePointProbeFromField {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GenerateSinglePointProbeFromField	
    }

    method move_location {} {
      set $this-moveto "location"
      $this-c needexecute
    }

    method move_node {} {
      set $this-moveto "node"
      $this-c needexecute
    }

    method move_elem {} {
      set $this-moveto "elem"
      $this-c needexecute
    }

    method move_center {} {
      set $this-moveto "center"
      $this-c needexecute
    }


    method changevalue {} {
      if { [set $this-show-value] } {
          $this-c needexecute
      } else {
          set $this-value ""
      }
    }

    method changenode {} {
      if { [set $this-show-node] } {
          $this-c needexecute
      } else {
          set $this-node ""
      }
    }

    method changeelem {} {
      if { [set $this-show-elem] } {
          $this-c needexecute
      } else {
          set $this-elem ""
      }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        build_ui $w

        makeSciButtonPanel $w $w $this -no_execute "\"Reset\" \"$this move_center\" \"\""
        moveToCursor $w
    }

    method build_ui { w } {
        global $this-main_frame
        set $this-main_frame $w
        
        sci_frame $w.f
        sci_frame $w.f.g
        sci_frame $w.f.g.labels
        sci_frame $w.f.g.entries
        sci_frame $w.f.g.entries.loc
        sci_frame $w.f.h

        sci_label $w.f.g.labels.location -text "Location" -just left
        sci_entry $w.f.g.entries.loc.locx -width 10 -textvariable $this-locx
        sci_entry $w.f.g.entries.loc.locy -width 10 -textvariable $this-locy
        sci_entry $w.f.g.entries.loc.locz -width 10 -textvariable $this-locz

        sci_checkbutton $w.f.g.labels.value -text "Value" -just left \
            -variable $this-show-value -command "$this changevalue"
        sci_entry $w.f.g.entries.value -width 40 -state disabled -textvariable $this-value

        sci_checkbutton  $w.f.g.labels.node -text "Node" -just left \
            -variable $this-show-node -command "$this changenode"
        sci_entry $w.f.g.entries.node -width 10 -textvariable $this-node

        sci_checkbutton $w.f.g.labels.elem -text "Elem" -just left \
            -variable $this-show-elem -command "$this changeelem"
        sci_entry $w.f.g.entries.elem -width 10 -textvariable $this-elem

        pack  $w.f.g.labels.location $w.f.g.labels.value \
            $w.f.g.labels.node $w.f.g.labels.elem\
          -side top -anchor w

        pack $w.f.g.entries.loc.locx $w.f.g.entries.loc.locy $w.f.g.entries.loc.locz \
          -side left -anchor n -expand yes -fill x

        pack $w.f.g.entries.loc -side top -expand yes -fill x
        pack $w.f.g.entries.value $w.f.g.entries.node $w.f.g.entries.elem \
          -side top -anchor w

        pack $w.f.g.labels $w.f.g.entries -side left


        sci_scale $w.f.slide -orient horizontal -label "GenerateSinglePointProbeFromField Size" -from 0 -to 100 -showvalue true \
             -variable $this-probe_scale -resolution 0.25 -tickinterval 25
        set $w.f.slide $this-probe_scale

        sci_label $w.f.h.label1 -text "Label" -just left
        sci_entry $w.f.h.entry1 -textvariable $this-label -width 20
        sci_frame $w.f.h.color

        addColorSelection $w.f.h.color "Color" $this-color "color_change"
        
        grid $w.f.h.label1 -row 0 -column 0 -sticky w
        grid $w.f.h.entry1 -row 0 -column 1 -sticky w
        grid $w.f.h.color -row 1 -column 0 -columnspan 2 -sticky w

        bind $w.f.slide <ButtonRelease> "$this-c needexecute"
        bind $w.f.slide <B1-Motion> "$this-c needexecute"

        pack $w.f.h $w.f.slide $w.f.g -side bottom -expand yes -fill x

        pack $w.f -side top -expand yes -fill both -padx 5 -pady 5

        bind $w.f.g.entries.loc.locx <KeyPress-Return> "$this move_location"
        bind $w.f.g.entries.loc.locy <KeyPress-Return> "$this move_location"
        bind $w.f.g.entries.loc.locz <KeyPress-Return> "$this move_location"
        bind $w.f.g.entries.node <KeyPress-Return> "$this move_node"
        bind $w.f.g.entries.elem <KeyPress-Return> "$this move_elem"
  
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




