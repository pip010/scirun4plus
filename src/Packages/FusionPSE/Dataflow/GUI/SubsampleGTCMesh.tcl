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


#    File   : SubsampleGTCMesh.tcl
#    Author : Allen Sanderson
#             SCI Institute
#             University of Utah
#    Date   : March 2008

# This GUI interface is for selecting values for sub sampling a
# topologically structured field.

catch {rename FusionPSE_Fields_SubsampleGTCMesh ""}

itcl::class FusionPSE_Fields_SubsampleGTCMesh {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name SubsampleGTCMesh
        set_defaults
    }

    method set_defaults {} {
	global power_app_command
	set    power_app_command ""

	for {set i 0} {$i < 2} {incr i 1} {
	    if { $i == 0 } {
		set index i
	    } elseif { $i == 1 } {
		set index j
	    }

	    global $this-dim-$index
	    trace variable $this-dim-$index w "$this update_set_size_callback"

	    global $this-start-$index
	    trace variable $this-start-$index w "$this update_SliderEntry4_callback"
	}
    }

    method set_power_app_cmd { cmd } {
	global power_app_command
	set power_app_command $cmd
    }

    method ui {} {

	set tmp 0.0

        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }

        toplevel $w

	frame $w.main

	frame $w.main.l
	label $w.main.l.direction -text "Index"  -width 5 -anchor w -just left
	label $w.main.l.start     -text "Start"  -width 5 -anchor w -just left
	label $w.main.l.stop      -text "Stop"   -width 5 -anchor w -just left

	pack $w.main.l.direction -side left -padx  20
	pack $w.main.l.start     -side left -padx  70
	pack $w.main.l.stop      -side left -padx 110

	for {set i 0} {$i < 2} {incr i 1} {
	    if { $i == 0 } {
		set index i
	    } elseif { $i == 1 } {
		set index j
	    }

	    global $this-dim-$index
	    global $this-start-$index
	    global $this-start2-$index
	    global $this-stop-$index
	    global $this-stop2-$index

	    set start_val 1
	    set stop_val [expr [set $this-dim-$index] - 2]

	    frame $w.main.$index

	    label $w.main.$index.l -text " $index :" \
		-width 3 -anchor w -just left

	    pack $w.main.$index.l -side left

	    scaleEntry4 $w.main.$index.start \
		0 $stop_val 200 \
		$this-start-$index $this-start2-$index $index

	    scaleEntry2 $w.main.$index.stop \
		$start_val [expr [set $this-dim-$index] - 1] 200 \
		$this-stop-$index $this-stop2-$index

	    pack $w.main.$index.l $w.main.$index.start $w.main.$index.stop \
		    -side left
	}
	
	pack $w.main.l $w.main.i $w.main.j -side top -padx 10 -pady 5

	pack $w.main -side top -fill x -expand 1
	
	global power_app_command

	if { [in_power_app] } {
	    makeSciButtonPanel $w $w $this -no_execute -no_close -no_find \
		"\"Close\" \"wm withdraw $w; $power_app_command\" \"Hides this GUI\""
	} else {
	    makeSciButtonPanel $w $w $this
	}
	 
	moveToCursor $w
    }

    method scaleEntry2 { win start stop length var1 var2 } {
	frame $win 
	pack $win -side top -padx 5

	scale $win.s -from $start -to $stop -length $length \
	    -variable $var1 -orient horizontal -showvalue false \
	    -command "$this updateSliderEntry $var1 $var2"

	entry $win.e -width 4 -text $var2

	bind $win.e <KeyRelease> "$this manualSliderEntry $start $stop $var1 $var2"

	pack $win.s -side left
	pack $win.e -side bottom -padx 5
    }

    method updateSliderEntry {var1 var2 someUknownVar} {
	set $var2 [set $var1]
    }

    method manualSliderEntry { start stop var1 var2 } {

	if { ![string is integer [set $var2]] } {
	    set $var2 [set $var1] }

	if { [set $var2] < $start } {
	    set $var2 $start }
	
	if { [set $var2] > $stop } {
	    set $var2 $stop }
	
	set $var1 [set $var2]
    }


    method scaleEntry4 { win start stop length var1 var2 index } {
	frame $win 
	pack $win -side top -padx 5

	scale $win.s -from $start -to $stop -length $length \
	    -variable $var1 -orient horizontal -showvalue false \
	    -command "$this updateSliderEntry4 $index"

	entry $win.e -width 4 -text $var2

	bind $win.e <KeyRelease> \
	    "$this manualSliderEntry4 $start $stop $var1 $var2 $index"

	pack $win.s -side left
	pack $win.e -side bottom -padx 5
    }

    method update_SliderEntry4_callback { name1 name2 op } {

	updateSliderEntry4 i 0
	updateSliderEntry4 j 0
    }

    method updateStopSlider { index } {

	global $this-start-$index
	global $this-start2-$index
	global $this-stop-$index
	global $this-stop2-$index
	global $this-dim-$index

	set w .ui[modname]

	if [ expr [winfo exists $w] ] {

	    # Update the sliders to have the new end values.

	    set start_val [expr [set $this-start-$index] + 1]
	    set stop_val  [expr [set $this-dim-$index] - 2]

	    $w.main.$index.start.s configure -from 0 -to $stop_val
	    $w.main.$index.stop.s configure \
		-from $start_val -to [expr [set $this-dim-$index] - 1]

	    bind $w.main.$index.start.e <KeyRelease> \
		"$this manualSliderEntry4 0 $stop_val $this-start-$index $this-start2-$index $index"
	    bind $w.main.$index.stop.e  <KeyRelease> \
		"$this manualSliderEntry $start_val [expr [set $this-dim-$index] - 1] $this-stop-$index $this-stop2-$index"
	}
    }

    method updateSliderEntry4 { index someUknownVar } {

	global $this-start-$index
	global $this-start2-$index
	global $this-stop-$index
	global $this-stop2-$index

	updateStopSlider $index

	set $this-start2-$index [set $this-start-$index]
	set $this-stop2-$index  [set $this-stop-$index]
    }

    method manualSliderEntry4 { start stop var1 var2 index } {

	if { ![string is integer [set $var2]] } {
	    set $var2 [set $var1] }

	if { [set $var2] < $start } {
	    set $var2 $start }
	
	if { [set $var2] > $stop } {
	    set $var2 $stop }
	
	set $var1 [set $var2]

	updateSliderEntry4 $index 0
    }

    method update_index { } {
	for {set i 0} {$i < 2} {incr i 1} {
	    if { $i == 0 } {
		set index i
	    } elseif { $i == 1 } {
		set index j
	    }

	    global $this-start-$index
	    global $this-start2-$index
	    global $this-stop-$index
	    global $this-stop2-$index

	    set $this-start2-$index  [set $this-start-$index]
	    set $this-stop2-$index   [set $this-stop-$index]
	}
    }

    method update_set_size_callback { name1 name2 op } {
	set_size
    }

    method set_size { } {
	global $this-i-dim
	global $this-j-dim

	set w .ui[modname]

	if [ expr [winfo exists $w] ] {
	    pack forget $w.main.i
	    pack forget $w.main.j
	    
	    pack $w.main.l $w.main.i $w.main.j -side top -padx 10 -pady 5
	}

	for {set i 0} {$i < 3} {incr i 1} {
	    if { $i == 0 } {
		set index i
	    } elseif { $i == 1 } {
		set index j
	    }

	    global $this-start-$index
	    global $this-start2-$index
	    global $this-stop-$index
	    global $this-stop2-$index

	    set stop_val1 [expr [set $this-dim-$index] - 1]
	    set stop_val2 [expr [set $this-dim-$index] - 2]

	    if [ expr [winfo exists $w] ] {

		# Update the sliders to the new bounds.
		$w.main.$index.start.s configure -from 0 -to $stop_val2
		$w.main.$index.stop.s  configure -from 0 -to $stop_val1

		bind $w.main.$index.start.e <KeyRelease> \
		    "$this manualSliderEntry4 0 $stop_val2 $this-start-$index $this-start2-$index $index"
		bind $w.main.$index.stop.e  <KeyRelease> \
		    "$this manualSliderEntry  1 $stop_val1 $this-stop-$index $this-stop2-$index"
	    }

	    # Update the stop values to be at the initials values.
	    set $this-start-$index 0	    
	    set $this-stop-$index  $stop_val1

	    # Update the text values.
	    set $this-start2-$index  [set $this-start-$index]
	    set $this-stop2-$index   [set $this-stop-$index]
	}
    }
}
