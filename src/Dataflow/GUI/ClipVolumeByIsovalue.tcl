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



catch {rename ClipVolumeByIsovalue ""}

itcl::class SCIRun_NewField_ClipVolumeByIsovalue {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name ClipVolumeByIsovalue
	set_defaults
    }

    method set_defaults {} {
	global $this-continuous
	set $this-continuous 0

	trace variable $this-update_type w "$this update_type_callback"
	trace variable $this-isoval-max w "$this update_minmax_callback"
    }
	
    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        sci_frame       $w.function
        sci_radiobutton $w.function.lte -text "Less Than"    -value 1 -variable $this-lte
        sci_radiobutton $w.function.gte -text "Greater Than" -value 0 -variable $this-lte

        pack $w.function.lte $w.function.gte -side left -padx 4

        pack $w.function -side top

        sci_frame $w.isoslider

        scaleEntry2 $w.isoslider.isoval \
            [set $this-isoval-min] [set $this-isoval-max] \
             4c $this-isoval $this-isoval-typed

        sci_labeledframe $w.isoslider.opt -labelpos nw -labeltext "Options"
        set opt [$w.isoslider.opt childsite]

        sci_optionmenu $opt.update -labeltext "Update:" \
          -labelpos w -command "$this set_update_type $opt.update"
        $opt.update insert end "On Release" Manual Auto
        $opt.update select [set $this-update_type]

        global $this-update
        set $this-update $opt.update

        pack $opt.update -side top -anchor w -pady 25

        pack $w.isoslider.isoval $w.isoslider.opt -side top -anchor w -fill x

        pack $w.isoslider -side top -pady 4

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }

    method update_type_callback { name1 name2 op } {
	set window .ui[modname]
	if {[winfo exists $window]} {
	    [set $this-update] select [set $this-update_type]
	}
    }

    method set_isoval {} {
	global $this-update

	set type [[set $this-update] get]

	if { $type == "On Release" } {
	    eval "$this-c needexecute"
	}
    }
    method set_update_type { w } {
	global $w
	global $this-continuous
	global $this-update_type

	set $this-update_type [$w get]
	if { [set $this-update_type] == "Auto" } {
	    set $this-continuous 1
	} else {
	    set $this-continuous 0
	}
    }

    method update_minmax_callback { name1 name2 op } {
	set_min_max
    }

    method set_min_max { } {
	set w .ui[modname]
	global $this-isoval-min
	global $this-isoval-max

	set min [set $this-isoval-min]
	set max [set $this-isoval-max]

	if [ expr [winfo exists $w] ] {

            if { $max == $min } {
                set max [expr $min + 1]
            }

	    set lg [expr floor( log10($max-$min) ) ]
	    set range [expr pow(10.0, $lg )]

	    set scale 1.0

	    if { $lg > 5.0 } {
		set scale [expr pow(10.0, $lg-5 )]
	    }

	    set win $w.isoslider.isoval

	    $win.l.s configure -from $min -to $max
	    $win.l.s configure -resolution [expr $range/(1.0e4*$scale)]
	    $win.l.s configure -tickinterval [expr ($max - $min)]

	    bind $win.r.e <Return> "$this manualSliderEntryReturn \
             $min $max $this-isoval $this-isoval-typed"
	    bind $win.r.e <KeyRelease> "$this manualSliderEntry \
             $min $max $this-isoval $this-isoval-typed"
	}
    }

    method scaleEntry2 { win start stop length var_slider var_typed } {
        sci_frame $win 

        sci_frame $win.l
        sci_frame $win.r
	
        if { $start == $stop } { 
            set stop [expr $start + 1]
        }

        set lg [expr floor( log10($stop-$start) ) ]
        set range [expr pow(10.0, $lg )]

        set scale 1.0

        if { $lg > 5.0 } {
            set scale [expr pow(10.0, $lg-5 )]
        }

        sci_scale $win.l.s \
          -from $start -to $stop \
          -length $length \
          -variable $var_slider -orient horizontal -showvalue false \
          -command "$this updateSliderEntry $var_slider $var_typed" \
          -resolution [expr $range/(1.0e4*$scale)] \
          -tickinterval [expr ($stop - $start)]

        sci_entry $win.r.e -width 7 -text $var_typed

        bind $win.l.s <ButtonRelease> "$this set_isoval"

        bind $win.r.e <Return> "$this manualSliderEntryReturn \
                   $start $stop $var_slider $var_typed"
        bind $win.r.e <KeyRelease> "$this manualSliderEntry \
                   $start $stop $var_slider $var_typed"

        pack $win.l.s -side top -expand 1 -fill x -padx 5
        pack $win.r.e -side top -padx 5 -pady 3
        pack $win.l -side left -expand 1 -fill x
        pack $win.r -side right -fill y
    }

    method updateSliderEntry {var_slider var_typed someUknownVar} {
      global $this-continuous
            global $var_slider
            global $var_typed
      set $var_typed [set $var_slider]

      if { [set $this-continuous] == 1.0 } {
          eval "$this-c needexecute"
      } elseif { [set $this-update_type] == "Auto" } {
          set $this-continuous 1
      }
    }

    method manualSliderEntryReturn { start stop var_slider var_typed } {
      # Since the user has typed in a value and hit return, we know
      # they are done and if their value is not valid or within range,
      # we can change it to be either the old value, or the min or max
      # depending on what is appropriate.
        if { ![string is double [set $var_typed]] } {
            set $var_typed [set $var_slider] 
        }

      if {[set $var_typed] < $start} {
          set $var_typed $start
      } elseif {[set $var_typed] > $stop} {
          set $var_typed $stop
      }

      # Force the update to be manual
      global $this-continuous
      set continuous [set $this-continuous]
      
      set $this-continuous 0
      
      set $var_slider [set $var_typed]
      
      set $this-continuous $continuous
    }
    
    method manualSliderEntry { start stop var_slider var_typed } {
	# Evaluate as the user types in an isoval but never change the value
	# they are typing in because they might not be done. Only update the
	# actual isoval when user has typed in a double and it is within range.
	
 	set var_new [set $var_slider]

 	# only update the value if it evaluates to a double 
	# and is within range
 	if {[string is double [set $var_typed]] && 
 	    $start <= [set $var_typed] && 
 	    [set $var_typed] <= $stop} {
 	    set var_new [set $var_typed]
 	}
	
	# Force the update to be manual
  	global $this-continuous
  	set continuous [set $this-continuous]
	
  	set $this-continuous 0
	
  	set $var_slider $var_new
	
  	set $this-continuous $continuous
    }
}
