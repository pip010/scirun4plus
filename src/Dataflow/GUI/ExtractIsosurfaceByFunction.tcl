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


#    File   : ExtractIsosurfaceByFunction.tcl
#    Author : Allen Sanderson
#             SCI Institute
#             University of Utah
#    Date   : March 2006

# This GUI interface is for selecting an axis and index for sub sampling a
# topologically structured field.

itcl::class SCIRun_NewField_ExtractIsosurfaceByFunction {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name ExtractIsosurfaceByFunction
        set_defaults
    }

    method set_defaults {} {
        global $this-continuous
        set $this-continuous 0

        global $this-active-slice-value-selection-tab
        set $this-active-slice-value-selection-tab 0
        set $this-help ""
    }

    method ui {} {
        global $this-function
        global $this-active-slice-value-selection-tab

        set oldmeth [set $this-active-slice-value-selection-tab]

        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }


        sci_toplevel $w
        sci_labeledframe $w.inf -labeltext "Extract Isosurface By Function"

        set infoframe [$w.inf childsite]
        sci_frame $infoframe.info
        pack $infoframe.info -side left
        set info $infoframe.info
        sci_label $info.info1 -text "Function: RESULT = function(DATA,POS,ELEMENT,INDEX,A,B,C,...)"
        sci_label $info.info2 -text "Input array: DATA (scalar/vector/tensor: data from field port) "
        sci_label $info.info3 -text "Input array: X, Y, Z (scalar: Cartensian coordinates of node/element)"
        sci_label $info.info4 -text "Input array: POS (vector: vector with node/element position)"
        sci_label $info.info6 -text "Input array: INDEX (scalar: number of the element)"
        sci_label $info.info7 -text "Input array: SIZE (scalar: number of elements)"
        sci_label $info.info8 -text "Input array: ELEMENT (element: object containing element properties)"

        grid $info.info1 -row 0 -column 0 -sticky w
        grid $info.info2 -row 1 -column 0 -sticky w
        grid $info.info3 -row 2 -column 0 -sticky w
        grid $info.info4 -row 3 -column 0 -sticky w
        grid $info.info6 -row 4 -column 0 -sticky w
        grid $info.info7 -row 5 -column 0 -sticky w
        grid $info.info8 -row 6 -column 0 -sticky w

        pack $w.inf -side top -anchor w -fill x


# Function definition
        sci_labeledframe $w.func -labelpos nw \
            -labeltext "Function Definition"
        set func [$w.func childsite]

        option add *textBackground white	
        sci_scrolledtext $func.text -height 60 -hscrollmode dynamic

        $func.text insert end [set $this-function]

        pack $w.func -side top -anchor w -expand true -fill both
        pack $func.text  -side top -expand true -fill both -padx 5


# Slice Value Selection Methods
        sci_labeledframe $w.slice -labelpos nw \
            -labeltext "Slice Value Selection"
        set isf [$w.slice childsite]

        global Color
        sci_tabnotebook $isf.tabs -raiseselect true -height 200 \
            -backdrop $Color(Basecolor)
        pack $isf.tabs -side top -fill x -expand 1


###### Slice Value using slider
        set sliceslider [$isf.tabs add -label "Slider" \
               -command "set $this-active-slice-value-selection-tab 0"]

        scaleEntry2 $sliceslider.sliceval \
            [set $this-slice-value-min] [set $this-slice-value-max] \
             4c $this-slice-value $this-slice-value-typed

        sci_labeledframe $sliceslider.opt -labelpos nw -labeltext "Options"
        set opt [$sliceslider.opt childsite]
        
        sci_optionmenu $opt.update -labeltext "Update:" \
          -labelpos w -command "$this set_update_type $opt.update"
        $opt.update insert end "On Release" Manual Auto

        $opt.update select [set $this-update_type]

        global $this-update
        set $this-update $opt.update

        pack $opt.update -side top -anchor w -pady 25

        pack $sliceslider.sliceval $sliceslider.opt -side top -anchor w -fill x

###### Slice Value using quantity	
        set slicequant [$isf.tabs add -label "Quantity" \
               -command "set $this-active-slice-value-selection-tab 1"]
        

###### Save the sliceval-quantity since the iwidget resets it
        global $this-slice-value-quantity
        set quantity [set $this-slice-value-quantity]
        sci_spinint $slicequant.q -labeltext "Number of evenly-spaced slices: " \
            -range {0 100} -step 1 \
            -textvariable $this-slice-value-quantity \
            -width 10 -fixed 10 -justify right
        
        $slicequant.q delete 0 end
        $slicequant.q insert 0 $quantity

        sci_frame $slicequant.f
        sci_label $slicequant.f.l -text "List of Slice Values:"
        sci_entry $slicequant.f.e -width 40 -text $this-quantity-list -state disabled
        pack $slicequant.f.l $slicequant.f.e -side left -fill both -expand true

        sci_frame $slicequant.m
        sci_radiobutton $slicequant.m.f -text "Field MinMax" \
          -variable $this-quantity-range -value "field" \
          -command "$this-c needexecute"
        sci_radiobutton $slicequant.m.m -text "Manual" \
          -variable $this-quantity-range -value "manual" \
          -command "$this-c needexecute"

        sci_frame $slicequant.m.t 
        sci_label $slicequant.m.t.minl -text "Min"
        sci_entry $slicequant.m.t.mine -width 6 -text $this-quantity-min
        sci_label $slicequant.m.t.maxl -text "Max"
        sci_entry $slicequant.m.t.maxe -width 6 -text $this-quantity-max
        bind $slicequant.m.t.mine <Return> "$this-c needexecute"
        bind $slicequant.m.t.maxe <Return> "$this-c needexecute"
        pack $slicequant.m.t.minl $slicequant.m.t.mine $slicequant.m.t.maxl $slicequant.m.t.maxe \
          -side left -fill x -expand true

        pack $slicequant.m.f -side top -anchor w
        pack $slicequant.m.m $slicequant.m.t -side left -anchor w

        sci_frame $slicequant.t
        sci_radiobutton $slicequant.t.e -text "Exclusive" \
          -variable $this-quantity-clusive -value "exclusive" \
          -command "$this-c needexecute"
        sci_radiobutton $slicequant.t.i -text "Inclusive" \
          -variable $this-quantity-clusive -value "inclusive" \
          -command "$this-c needexecute"

        pack $slicequant.t.e $slicequant.t.i -side left -anchor w

        pack $slicequant.q $slicequant.m $slicequant.t -side top -expand true -fill x -pady 5

        pack $slicequant.f -fill x

###### Slice Value using list
        set slicelist [$isf.tabs add -label "List" \
             -command "set $this-active-slice-value-selection-tab 2"]

        
        sci_frame $slicelist.f
        sci_label $slicelist.f.l -text "List of Slice Values:"
        sci_entry $slicelist.f.e -width 40 -text $this-slice-value-list
        bind $slicelist.f.e <Return> "$this-c needexecute"
        pack $slicelist.f.l $slicelist.f.e -side left -fill both -expand true
        pack $slicelist.f -fill x


###### Slice Value using matrix
        set slicematrix [$isf.tabs add -label "Matrix" \
               -command "set $this-active-slice-value-selection-tab 3"]

        sci_frame $slicematrix.f
        sci_label $slicematrix.f.l -text "List of Slice Values:"
        sci_entry $slicematrix.f.e -width 40 -text $this-matrix-list -state disabled
        pack $slicematrix.f.l $slicematrix.f.e -side left -fill both -expand true
        pack $slicematrix.f -fill x


# Pack the Slice Value Selection Tabs

        $isf.tabs view $oldmeth
        $isf.tabs configure -tabpos "n"

        pack $isf.tabs -side top

        pack $w.slice -side top -anchor w -expand true -fill x
  
        sci_button $w.help -text "Available Functions" -command "$this showhelp"
        pack $w.help -side top -anchor e

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
    
    method showhelp { } {
      # Create a unique name for the file selection window
      set w [format "%s-functionhelp" .ui[modname]]

      if { [winfo exists $w] } {
        if { [winfo ismapped $w] == 1} {
          raise $w
        } else {
          wm deiconify $w
        }
	    	return
      }
	
      sci_toplevel $w -class TkFDialog
            
      sci_labeledframe $w.hf -labeltext "Parser Help"
      set help [$w.hf childsite]
      option add *textBackground white	
      sci_scrolledhtml $help.help -height 60 -hscrollmode dynamic -width 500p -height 300p        
      
      set helpfile [file join [netedit getenv SCIRUN_SRCDIR] Dataflow GUI ArrayMathFunctionHelp.html]
      $help.help import $helpfile
      
      pack $help.help -side top -anchor w -fill both -expand yes
      pack $w.hf -side top -anchor w -fill both -expand yes
    } 

    method update_text {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
          set func [$w.func childsite]
          set $this-function [$func.text get 1.0 end]
        }
    }

    method set-slice-value {} {
	global $this-update

	set type [[set $this-update] get]
	if { $type == "On Release" } {
	    eval "$this-c needexecute"
	}
    }
    
    method set-slice-quant-list { vals } {
	global $this-quantity-list
	
	set $this-quantity-list $vals
    }
    
    method set-slice-matrix-list { vals } {
	global $this-matrix-list
	
	set $this-matrix-list $vals
    }

    method update_type_callback { name1 name2 op } {
        set tmp [set $this-update_type]
        if { $tmp == "on release" } { set $this-update_type "On Release" }
	set window .ui[modname]
	if {[winfo exists $window]} {
	    set opt [$window.f.opt childsite]
	    $opt.update select [set $this-update_type]
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
	global $this-slice-value-min
	global $this-slice-value-max

	set min [set $this-slice-value-min]
	set max [set $this-slice-value-max]

	if [ expr [winfo exists $w] ] {
	    set lg [expr floor( log10($max-$min) ) ]
	    set range [expr pow(10.0, $lg )]

	    set scale 1.0

	    if { $lg > 5.0 } {
		set scale [expr pow(10.0, $lg-5 )]
	    }

	    set win $w.slice.childsite.tabs.canvas.notebook.cs.page1.cs.sliceval

	    $win.l.s configure -from $min -to $max
	    $win.l.s configure -resolution [expr $range/(1.0e4*$scale)]
	    $win.l.s configure -tickinterval [expr ($max - $min)]

	    bind $win.r.e <Return> "$this manualSliderEntryReturn \
             $min $max $this-slice-value $this-slice-value-typed"
	    bind $win.r.e <KeyRelease> "$this manualSliderEntry \
             $min $max $this-slice-value $this-slice-value-typed"
	}
    }

    method scaleEntry2 { win start stop length var_slider var_typed } {
	sci_frame $win 

	sci_frame $win.l
	sci_frame $win.r
	
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

	bind $win.l.s <ButtonRelease> "$this set-slice-value"

	bind $win.r.e <Return> "$this manualSliderEntryReturn \
             $start $stop $var_slider $var_typed"
	bind $win.r.e <KeyRelease> "$this manualSliderEntry \
             $start $stop $var_slider $var_typed"

	pack $win.l.s -side top -expand true -fill x -padx 5
	pack $win.r.e -side top -padx 5 -pady 3
	pack $win.l -side left -expand true -fill x
	pack $win.r -side right -fill y
    }

    method updateSliderEntry {var_slider var_typed someUknownVar} {
	global $this-continuous
	global $this-update_type
        global $var_typed
	set $var_typed [set $var_slider]
	
	if { [set $this-continuous] == 1.0 } {
	    eval "$this-c needexecute"
	} elseif { [set $this-update_type] == "Auto" } {
	    set $this-continuous 1
	}
    }

    method manualSliderEntryReturn { start stop var_slider var_typed } {
	# Because the user has typed in a value and hit return, we know
	# they are done and if their value is not valid or within range,
	# we can change it to be either the old value, or the min or max
	# depending on what is appropriate.
  	if { ![string is integer [set $var_typed]] } {
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

	if { [set $this-update_type] == "On Release" ||
	     [set $this-update_type] == "Auto" } {
	    eval "$this-c needexecute"
	}
    }


    method manualSliderEntry { start stop var_slider var_typed } {
	# Evaluate as the user types in an sliceval but never change the value
	# they are typing in because they might not be done. Only update the
	# actual sliceval when user has typed in a double and it is within range.
	
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

