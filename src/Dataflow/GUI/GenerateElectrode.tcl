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


itcl::class SCIRun_NewField_GenerateElectrode {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GenerateElectrode	
    }
		
		method move_num_points {} {
      set $this-moveto "num_points"
      $this-c needexecute
    }

		
		method move_default {} {
      set $this-moveto "default"
      $this-c needexecute
    }
		
		method add_point {} {
      set $this-moveto "add_point"
      $this-c needexecute
    }
		
		method remove_point {} {
      set $this-moveto "remove_point"
      $this-c needexecute
    }
		
		method enable_options {} {
			
		 set w .ui[modname]
		 
		 $w.f.g.labels.width configure -foreground black
		 $w.f.g.entries.width configure -state normal -foreground black
		 $w.f.proj.l configure	-foreground black
		 $w.f.proj.f1.b configure -state normal -foreground black
		 $w.f.proj.f2.b configure -state normal -foreground black
		 $w.f.proj.f3.b configure -state normal -foreground black
	 }
	 
	 method disable_options {} {
	 
		 set w .ui[modname]		 
		 set tcolor "#505050"
	 
		 $w.f.g.labels.width configure -foreground $tcolor
		 $w.f.g.entries.width configure -state disabled -foreground $tcolor
		 $w.f.proj.l configure	-foreground $tcolor
		 $w.f.proj.f1.b configure -state disabled -foreground $tcolor
		 $w.f.proj.f2.b configure -state disabled -foreground $tcolor
		 $w.f.proj.f3.b configure -state disabled -foreground $tcolor
	  }
		
		method wire_type {} {
	 
		 global $this-electrode_type		 
		 set $this-electrode_type "wire"
		 		 
		 disable_options	 
		 
		 $this-c needexecute
	 
	 }
	 
	 method planar_type {} {
	 
		 global $this-electrode_type		 
		 set $this-electrode_type "planar"
		 
		 enable_options
		 
		 $this-c needexecute
	 
	 }
	 
	 method check_status {} {
	 
		 set type [set $this-electrode_type]
		 
		 if {$type=="wire"} {disable_options}
		 if {$type=="planar"} {enable_options}		 
	 }

		
		method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        build_ui $w

        makeSciButtonPanel $w $w $this "\"Reset\" \"$this move_default\" \"\""
        moveToCursor $w
    }

    method build_ui { w } {
        global $this-main_frame
        set $this-main_frame $w
				
        sci_frame $w.f
				sci_frame $w.f.b
        sci_frame $w.f.g
				sci_frame $w.f.g.labels
        sci_frame $w.f.g.entries
        sci_frame $w.f.h
				sci_frame $w.f.p
				sci_frame $w.f.slide
				
				sci_frame $w.type -relief groove -border 2
        sci_label $w.type.l -text "Electrode Type"
        pack $w.type -side top -padx 2 -pady 2 -fill both
        pack $w.type.l -side top -fill x
			
        sci_frame $w.type.f1 -relief flat
        pack $w.type.f1 -side top -expand yes -fill x
        sci_radiobutton $w.type.f1.b -text "Wire Electrode" \
            -variable $this-electrode_type -value "wire" \
						-command "$this wire_type"
        pack $w.type.f1.b -side left

        sci_frame $w.type.f2 -relief flat
        pack $w.type.f2 -side top -expand yes -fill x
        sci_radiobutton $w.type.f2.b -text "Planar Electrode" \
          -variable $this-electrode_type -value "planar" \
					-command "$this planar_type"
        pack $w.type.f2.b -side left
				
				sci_checkbutton  $w.f.b.use_field -text "Use Field Nodes" -just left \
						-onvalue 1 -offvalue 0 -variable $this-use-field -command "$this-c needexecute"
				sci_checkbutton  $w.f.b.move_all -text "Move All Nodes" -just left \
						-onvalue 1 -offvalue 0 -variable $this-move-all -command "$this-c needexecute"
												
				grid $w.f.b.use_field -row 0 -column 0 -columnspan 2 -sticky w
				grid $w.f.b.move_all -row 1 -column 0 -columnspan 2 -sticky w
				pack $w.f.b -side top			
							
				
				sci_label $w.f.g.labels.length -text "Length of Electrode" -just left
				sci_entry $w.f.g.entries.length -width 10 -textvariable $this-length
				pack  $w.f.g.labels.length -side top -anchor w
				pack $w.f.g.entries.length -side left -anchor n -expand yes -fill x
				
				
				grid $w.f.g.labels.length -row 1 -column 0 -sticky w
        grid $w.f.g.entries.length -row 1 -column 1 -sticky w
				
				sci_label $w.f.g.labels.width -text "Width of Electrode" -just left 
				$w.f.g.labels.width configure
				sci_entry $w.f.g.entries.width -width 10 -textvariable $this-width
				pack  $w.f.g.labels.width -side top -anchor w
				pack $w.f.g.entries.width -side left -anchor n -expand yes -fill x
				
				grid $w.f.g.labels.width -row 2 -column 0 -sticky w
        grid $w.f.g.entries.width -row 2 -column 1 -sticky w
				
				sci_label $w.f.g.labels.thick -text "Thickness of Electrode" -just left
				sci_entry $w.f.g.entries.thick -width 10 -textvariable $this-thick
				pack  $w.f.g.labels.thick -side top -anchor w
				pack $w.f.g.entries.thick -side left -anchor n -expand yes -fill x
				
				grid $w.f.g.labels.thick -row 3 -column 0 -sticky w
        grid $w.f.g.entries.thick -row 3 -column 1 -sticky w
				
				sci_label $w.f.g.labels.num_points -text "Number of Control Points" -just left
				sci_entry $w.f.g.entries.num_points -width 10 -textvariable $this-num_points
				pack  $w.f.g.labels.num_points -side top -anchor w
				pack $w.f.g.entries.num_points -side left -anchor n -expand yes -fill x
				
				grid $w.f.g.labels.num_points -row 4 -column 0 -sticky w
        grid $w.f.g.entries.num_points -row 4 -column 1 -sticky w
				
				pack $w.f.g.labels $w.f.g.entries -side left
				
				sci_frame $w.f.proj -relief groove -border 2 
        sci_label $w.f.proj.l -text "Projection of Electrode Along Normal Vector"
				$w.f.proj.l configure
        pack $w.f.proj -side top -padx 2 -pady 2 -fill both
        pack $w.f.proj.l -side top -fill x
				
				
				sci_frame $w.f.proj.f1 -relief flat
        pack $w.f.proj.f1 -side top -expand yes -fill x
        sci_radiobutton $w.f.proj.f1.b -text "Positive" \
            -variable $this-project -value "positive" \
						-command "$this-c needexecute"
        pack $w.f.proj.f1.b -side left
				
				sci_frame $w.f.proj.f2 -relief flat
				pack $w.f.proj.f2 -side top -expand yes -fill x
        sci_radiobutton $w.f.proj.f2.b -text "Midway" \
            -variable $this-project -value "midway" \
						-command "$this-c needexecute"
        pack $w.f.proj.f2.b -side left
				
				sci_frame $w.f.proj.f3 -relief flat						
        pack $w.f.proj.f3 -side top -expand yes -fill x
        sci_radiobutton $w.f.proj.f3.b -text "Negative" \
            -variable $this-project -value "negative" \
						-command "$this-c needexecute"
        pack $w.f.proj.f3.b -side left
				
				pack $w.f.p -side bottom
				
				check_status
					
				sci_scale $w.f.slide.size -orient horizontal -label "Size of Point Widgets" -from 0 -to 20 -showvalue true \
             -variable $this-probe_scale -resolution 0.25 -tickinterval 5
        set $w.f.slide.size $this-probe_scale
				
				sci_scale $w.f.slide.res -orient horizontal -label "Resolution of the Electrode" -from 1 -to 50 -showvalue true \
             -variable $this-wire_res -resolution 1 -tickinterval 49
        set $w.f.slide.res $this-wire_res
      
        sci_frame $w.f.h.color

        addColorSelection $w.f.h.color "Color" $this-color "color_change"
				
				sci_button $w.f.p.add_point -text "Add Point" -command "$this add_point"
        #grid   $w.f.p.add_point -row 0 -column 0
								
				sci_button $w.f.p.remove_point -text "Remove Point" -command "$this remove_point"
        #grid   $w.f.p.remove_point -row 0 -column 1

        pack $w.f.p.add_point $w.f.p.remove_point -side top -anchor e
				
        
				
        grid $w.f.h.color -row 1 -column 0 -columnspan 2 -sticky w

        bind $w.f.slide.size <ButtonRelease> "$this-c needexecute"
				bind $w.f.slide.res <ButtonRelease> "$this-c needexecute"
        
        pack $w.f.h $w.f.slide $w.f.slide.size $w.f.slide.res $w.f.g -side bottom -expand yes -fill x

        pack $w.f -side top -expand yes -fill both -padx 5 -pady 5
				
				
				
				
				bind $w.f.g.entries.length <KeyPress-Return> "$this-c needexecute"
				bind $w.f.g.entries.num_points <KeyPress-Return> "$this move_num_points"
				bind $w.f.g.entries.width <KeyPress-Return> "$this-c needexecute"

        
  
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




