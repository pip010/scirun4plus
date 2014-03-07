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



catch {rename ExtractIsosurface ""}

itcl::class SCIRun_Visualization_ExtractIsosurface {
    inherit Module
    
     constructor { {args ""} } {
        eval configure $args
	set name ExtractIsosurface
	set_defaults
    }
    
    method set_defaults {} {
	global $this-continuous
	set $this-continuous 0

	trace variable $this-active_tab w "$this switch_to_active_tab"
	trace variable $this-update_type w "$this update_type_callback"
	trace variable $this-isoval-max w "$this update_minmax_callback"

	# SAGE vars
	global $this-visibility $this-value $this-scan
	global $this-bbox
	global $this-cutoff_depth 
	global $this-reduce
	global $this-all
	global $this-rebuild
	global $this-min_size
	global $this-poll

	set $this-visiblilty 0
	set $this-value 1
	set $this-scan 1
	set $this-bbox 1
	set $this-reduce 1
	set $this-all 0
	set $this-rebuild 0
	set $this-min_size 1
	set $this-poll 0

    }

    method raiseColor {col color colMsg} {
	 global $color
	 set window .ui[modname]
	 if {[winfo exists $window.color]} {
	     SciRaise $window.color
	     return
	 } else {
	     # makeColorPicker now creates the $window.color toplevel.
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

    method switch_to_active_tab {name1 name2 op} {
      set window .ui[modname]
      if {[winfo exists $window]} {
          set mf [$window.f.meth childsite]
          $mf.tabs view [set $this-active_tab]
      }
    }

    method ui {} {
      set w .ui[modname]
      if {[winfo exists $w]} {
          return
      }
      
      sci_toplevel $w
      sci_frame $w.f 
      pack $w.f -padx 2 -pady 2 -expand 1 -fill x
      set n "$this-c needexecute"

      set oldmeth [set $this-active-isoval-selection-tab]

    # Iso Value Selection Methods
      sci_labeledframe $w.f.iso -labelpos nw \
          -labeltext "Isovalue Selection"
      set isf [$w.f.iso childsite]
      global Color
      sci_tabnotebook $isf.tabs -raiseselect true -height 200 \
          -backdrop $Color(Basecolor)
      pack $isf.tabs -side top -fill x -expand 1
      pack $w.f.iso -side top -fill x -expand 1


    ###### Iso Value using slider
      set isoslider [$isf.tabs add -label "Slider" \
             -command "set $this-active-isoval-selection-tab 0"]

      scaleEntry2 $isoslider.isoval \
          [set $this-isoval-min] [set $this-isoval-max] \
           4c $this-isoval $this-isoval-typed

      sci_labeledframe $isoslider.opt -labelpos nw -labeltext "Options"
      set opt [$isoslider.opt childsite]

      sci_optionmenu $opt.update -labeltext "Update:" \
        -labelpos w -command "$this set_update_type $opt.update"
      $opt.update insert end "On Release" Manual Auto
      $opt.update select [set $this-update_type]

      global $this-update
      set $this-update $opt.update

      pack $opt.update -side top -anchor w -pady 25

      pack $isoslider.isoval $isoslider.opt -side top -anchor w -fill x


###### Iso Value using quantity	
      set isoquant [$isf.tabs add -label "Quantity" \
             -command "set $this-active-isoval-selection-tab 1"]
      

###### Save the isoval-quantity since the iwidget resets it
      global $this-isoval-quantity
      set quantity [set $this-isoval-quantity]
      sci_spinint $isoquant.q -labeltext "Number of evenly-spaced isovals: " \
          -range {0 100} -step 1 \
          -textvariable $this-isoval-quantity \
          -width 10 -fixed 10 -justify right
      
      $isoquant.q delete 0 end
      $isoquant.q insert 0 $quantity

      sci_frame $isoquant.f
      sci_label $isoquant.f.l -text "List of Isovals:"
      sci_entry $isoquant.f.e -width 40 -text $this-quantity-list -state disabled
      pack $isoquant.f.l $isoquant.f.e -side left -fill both -expand 1

      sci_frame $isoquant.m
      sci_radiobutton $isoquant.m.c -text "ColorMap MinMax" \
        -variable $this-quantity-range -value "colormap"
#        -command "$this-c needexecute"
      sci_radiobutton $isoquant.m.f -text "Field MinMax" \
        -variable $this-quantity-range -value "field"
#        -command "$this-c needexecute"
      sci_radiobutton $isoquant.m.m -text "Manual" \
        -variable $this-quantity-range -value "manual"
#        -command "$this-c needexecute"

      sci_frame $isoquant.m.t 
      sci_label $isoquant.m.t.minl -text "Min"
      sci_entry $isoquant.m.t.mine -width 6 -text $this-quantity-min
      sci_label $isoquant.m.t.maxl -text "Max"
      sci_entry $isoquant.m.t.maxe -width 6 -text $this-quantity-max
#      bind $isoquant.m.t.mine <Return> "$this-c needexecute"
#      bind $isoquant.m.t.maxe <Return> "$this-c needexecute"
      pack $isoquant.m.t.minl $isoquant.m.t.mine $isoquant.m.t.maxl $isoquant.m.t.maxe \
        -side left -fill x -expand 1

      pack $isoquant.m.c $isoquant.m.f -side top -anchor w
      pack $isoquant.m.m $isoquant.m.t -side left -anchor w

      sci_frame $isoquant.t
      sci_radiobutton $isoquant.t.e -text "Exclusive" \
        -variable $this-quantity-clusive -value "exclusive" 
#        -command "$this-c needexecute"
      sci_radiobutton $isoquant.t.i -text "Inclusive" \
        -variable $this-quantity-clusive -value "inclusive" 
#        -command "$this-c needexecute"

      pack $isoquant.t.e $isoquant.t.i -side left -anchor w

      pack $isoquant.q $isoquant.m $isoquant.t -side top -expand 1 -fill x -pady 5

      pack $isoquant.f -fill x


      ###### Iso Value using list
      set isolist [$isf.tabs add -label "List" \
           -command "set $this-active-isoval-selection-tab 2"]

      
      sci_frame $isolist.f
      sci_label $isolist.f.l -text "List of Isovals:"
      sci_entry $isolist.f.e -width 40 -text $this-isoval-list
#      bind $isolist.f.e <Return> "$this-c needexecute"
      pack $isolist.f.l $isolist.f.e -side left -fill both -expand 1
      pack $isolist.f -fill x


      ###### Iso Value using matrix
      set isomatrix [$isf.tabs add -label "Matrix" \
             -command "set $this-active-isoval-selection-tab 3"]

      sci_frame $isomatrix.f
      sci_label $isomatrix.f.l -text "List of Isovals:"
      sci_entry $isomatrix.f.e -width 40 -text $this-matrix-list -state disabled
      pack $isomatrix.f.l $isomatrix.f.e -side left -fill both -expand 1
      pack $isomatrix.f -fill x


      # Pack the Iso Value Selection Tabs

      $isf.tabs view $oldmeth
      $isf.tabs configure -tabpos "n"

      pack $isf.tabs -side top
      pack $w.f.iso -side top


      #  Options

      sci_labeledframe $w.f.options -labelpos nw -labeltext "Options"
      set opts [$w.f.options childsite]
      
      global $this-build_trisurf
      sci_checkbutton $opts.buildsurf -text "Build Output Field" \
        -variable $this-build_trisurf

      global $this-build_geom
      sci_frame $opts.buildgeom

      sci_checkbutton $opts.buildgeom.check -text "Build Output Geometry" \
        -variable $this-build_geom

      sci_checkbutton $opts.buildgeom.tranps -text "Enable Transparency (Geometry Only)" \
        -variable $this-transparency

      pack $opts.buildgeom.check $opts.buildgeom.tranps -side top -anchor w

      pack $opts.buildsurf $opts.buildgeom \
          -side top -anchor w

      addColorSelection $opts "Default Color" \
          $this-color "default_color_change"

      pack $w.f.options -side top -fill x -expand 1

      #  Methods
      sci_labeledframe $w.f.meth -labelpos nw \
          -labeltext "Computation Method"
      set mf [$w.f.meth childsite]

      sci_frame $mf.mc
      sci_radiobutton $mf.mc.r1 -text "Marching Cubes (Multi Threaded)" \
          -variable $this-num-threads -value 0 -command "$this select-alg"
      sci_radiobutton $mf.mc.r2 -text "Marching Cubes (Single Threaded)" \
          -variable $this-num-threads -value 1 -command "$this select-alg"

      pack $mf.mc.r1 $mf.mc.r2 -side top

      pack $mf.mc -side top -anchor w -expand y -fill x
      pack $w.f.meth -side top -fill x -expand 1

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }

    method set-isoval {} {
      global $this-update

      set type [[set $this-update] get]

      if { $type == "On Release" } {
          eval "$this-c needexecute"
      }
        }
        
    method set-isoquant-list { vals } {
      global $this-quantity-list
      
      set $this-quantity-list $vals
        }
        
    method set-isomatrix-list { vals } {
      global $this-matrix-list
      
      set $this-matrix-list $vals
        }
        
    method orient { tab page { val 4 }} {
      global $page
      global $tab
      
      $tab.tabs configure -tabpos [$page.orient get]
    }

    method select-alg {} {
      global $this-update

      set type [[set $this-update] get]

      if { $type != "Manual" } {
          eval "$this-c needexecute"
      }
        }

    method update_type_callback { name1 name2 op } {
      set window .ui[modname]
      if {[winfo exists $window]} {
          [set $this-update] select [set $this-update_type]
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

	    set win $w.f.iso.childsite.tabs.canvas.notebook.cs.page1.cs.isoval

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

        bind $win.l.s <ButtonRelease> "$this set-isoval"

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
