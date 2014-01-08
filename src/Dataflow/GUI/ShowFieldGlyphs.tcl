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


itcl::class SCIRun_Visualization_ShowFieldGlyphs {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
	set name ShowFieldGlyphs
	set_defaults
    }

    method set_defaults {} {

	trace variable $this-scalars_scaleNV w "$this new_scalars_scale"
	trace variable $this-scalars_thresholdNV w "$this new_scalars_threshold"
	trace variable $this-vectors_scaleNV w "$this new_vectors_scale"
	trace variable $this-vectors_thresholdNV w "$this new_vectors_threshold"
	trace variable $this-tensors_scaleNV w "$this new_tensors_scale"
	trace variable $this-tensors_thresholdNV w "$this new_tensors_threshold"
	trace variable $this-secondary_scaleNV w "$this new_secondary_scale"
	trace variable $this-tertiary_scaleNV w "$this new_tertiary_scale"

	trace variable $this-active_tab w "$this switch_to_active_tab"
	trace variable $this-scalars_has_data w "$this scalar_tab_changed"
	trace variable $this-vectors_has_data w "$this vector_tab_changed"
	trace variable $this-tensors_has_data w "$this tensor_tab_changed"
	trace variable $this-secondary_has_data w "$this secondary_tab_changed"
	trace variable $this-tertiary_has_data w "$this tertiary_tab_changed"

	# no C side component for these variables
	global $this-ss_slider
	set $this-ss_slider "not.set.yet"

	global $this-vs_slider
	set $this-vs_slider "not.set.yet"

	global $this-ts_slider
	set $this-ts_slider "not.set.yet"

	global $this-st_slider
	set $this-st_slider "not.set.yet"

	global $this-vt_slider
	set $this-vt_slider "not.set.yet"

	global $this-tt_slider
	set $this-tt_slider "not.set.yet"

	global $this-second_slider
	set $this-second_slider "not.set.yet"

	global $this-third_slider
	set $this-second_third "not.set.yet"
    }

    method new_scalars_scale {a1 a2 a3} {
      global $this-ss_slider
      set val [set $this-scalars_scaleNV]
      if {$val != -0.0} {
        upvar $this-ss_slider ss_slider
        if {[info exists $this-ss_slider] && [winfo exists $ss_slider]} {
            $ss_slider newvalue $val
        } else {
            set $this-scalars_scale [set $this-scalars_scaleNV]
        }
        set $this-scalars_scaleNV -0.0
      }
    }

    method new_scalars_threshold {a1 a2 a3} {
      global $this-st_slider
      set val [set $this-scalars_thresholdNV]
      if {$val != -0.0} {
        upvar $this-st_slider st_slider
        if {[info exists $this-st_slider] && [winfo exists $st_slider]} {
            $st_slider newvalue $val
        } else {
            set $this-scalars_threshold [set $this-scalars_thresholdNV]
        }
        set $this-scalars_thresholdNV -0.0
      }
    }
    
    method new_vectors_scale {a1 a2 a3} {
      global $this-vs_slider
      set val [set $this-vectors_scaleNV]
      if {$val != -0.0} {
        upvar $this-vs_slider vs_slider
        if {[info exists $this-vs_slider] && [winfo exists $vs_slider]} {
            $vs_slider newvalue $val
        } else {
            set $this-vectors_scale [set $this-vectors_scaleNV]
        }
        set $this-vectors_scaleNV -0.0
      }
    }


    method new_vectors_threshold {a1 a2 a3} {
      global $this-st_slider
      set val [set $this-vectors_thresholdNV]
      if {$val != -0.0} {
        upvar $this-st_slider st_slider
        if {[info exists $this-st_slider] && [winfo exists $st_slider]} {
            $st_slider newvalue $val
        } else {
            set $this-vectors_threshold [set $this-vectors_thresholdNV]
        }
        set $this-vectors_thresholdNV -0.0
      }
    }
    
    method new_tensors_scale {a1 a2 a3} {
      global $this-ts_slider
      set val [set $this-tensors_scaleNV]
      if {$val != -0.0} {
        upvar $this-ts_slider ts_slider
        if {[info exists $this-ts_slider] && [winfo exists $ts_slider]} {
            $ts_slider newvalue $val
        } else {
            set $this-tensors_scale [set $this-tensors_scaleNV]
        }
        set $this-tensors_scaleNV -0.0
      }
    }


    method new_tensors_threshold {a1 a2 a3} {
      global $this-st_slider
      set val [set $this-tensors_thresholdNV]
      if {$val != -0.0} {
        upvar $this-st_slider st_slider
        if {[info exists $this-st_slider] && [winfo exists $st_slider]} {
            $st_slider newvalue $val
        } else {
            set $this-tensors_threshold [set $this-tensors_thresholdNV]
        }
        set $this-tensors_thresholdNV -0.0
      }
    }
    method new_secondary_scale {a1 a2 a3} {
      global $this-second_slider
      set val [set $this-secondary_scaleNV]
      if {$val != -0.0} {
        upvar $this-second_slider second_slider
        if {[info exists $this-second_slider] && [winfo exists $second_slider]} {
            $second_slider newvalue $val
        } else {
            set $this-secondary_scale [set $this-secondary_scaleNV]
        }
        set $this-secondary_scaleNV -0.0
      }
    }
    
    method new_tertiary_scale {a1 a2 a3} {
      global $this-third_slider
      set val [set $this-tertiary_scaleNV]
      if {$val != -0.0} {
        upvar $this-third_slider third_slider
        if {[info exists $this-third_slider] && [winfo exists $third_slider]} {
            $third_slider newvalue $val
        } else {
            set $this-tertiary_scale [set $this-tertiary_scaleNV]
        }
        set $this-tertiary_scaleNV -0.0
      }
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


    method set_active_tab {act} {
	global $this-active_tab
	#puts stdout $act
	set $this-active_tab $act
    }

    method switch_to_active_tab {name1 name2 op} {
	#puts stdout "switching"
	set window .ui[modname]
	if {[winfo exists $window]} {
	    set dof [$window.options.disp.frame_title childsite]
	    $dof.tabs view [set $this-active_tab]
	}
    }

    # Text Tab
    method add_text_tab {dof} {
	set text [$dof.tabs add -label "Text" \
		-command "$this set_active_tab \"Text\""]
	sci_checkbutton $text.show_text \
		-text "Show Text" \
		-command "$this-c toggle_display_text" \
		-variable $this-text_on

	sci_frame $text.def_col -borderwidth 2

	sci_checkbutton $text.backfacecull \
	    -text "Cull backfacing text if possible" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_backface_cull

	sci_checkbutton $text.alwaysvisible \
	    -text "Text always visible (not hidden by faces)" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_always_visible

	sci_checkbutton $text.locations \
	    -text "Render indices as locations" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_render_locations

	sci_frame $text.show 
	sci_checkbutton $text.show.data \
	    -text "Show data values" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_show_data
	sci_checkbutton $text.show.nodes \
	    -text "Show node indices" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_show_nodes
	sci_checkbutton $text.show.edges \
	    -text "Show edge indices" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_show_edges
	sci_checkbutton $text.show.faces \
	    -text "Show face indices" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_show_faces
	sci_checkbutton $text.show.cells \
	    -text "Show cell indices" \
	    -command "$this-c rerender_text" \
	    -variable $this-text_show_cells

	make_labeled_radio $text.size \
	    "Text Size:" "$this-c rerender_text" left 5 \
	    $this-text_fontsize \
	    {{"XS" 0} {"S" 1} {"M" 2} {"L" 3} {"XL" 4}}

	pack $text.show.data $text.show.nodes $text.show.edges \
	    $text.show.faces $text.show.cells \
	    -side top -fill y -anchor w
	
	global $this-text_color_type
	make_labeled_radio $text.color \
	    "Text Coloring" "$this-c rerender_text" top 3 \
	    $this-text_color_type \
	    { {Text 0} {"Colormap Lookup" 1} {"RGB Conversion" 2} }

	addColorSelection $text.def_col "Text Color" $this-text_color \
	    "text_color_change"

	sci_frame $text.precision
	sci_label $text.precision.label -text "Text Precision  "
	sci_scale $text.precision.scale -orient horizontal \
	    -variable $this-text_precision -from 1 -to 16 \
	    -showvalue true -resolution 1
	bind $text.precision.scale <ButtonRelease> "$this-c rerender_text"
	pack $text.precision.label -side left -anchor s -fill x
	pack $text.precision.scale -side left -anchor n -fill x
	
	pack $text.show_text $text.show $text.backfacecull \
	    $text.alwaysvisible $text.locations \
	    $text.color $text.def_col $text.size \
	    $text.precision -side top -pady 2 -fill y -anchor w
    }


    # Scalar Tab
    method add_scalar_tab {dof} {

	set scalar [$dof.tabs add -label "Scalars" \
		-command "$this set_active_tab \"Scalars\""]
	sci_checkbutton $scalar.show_scalars \
		-text "Show Scalars" \
		-command "$this-c toggle_display_scalars" \
		-variable $this-scalars_on

	sci_checkbutton $scalar.transparency \
		-text "Enable Transparency" \
		-command "$this-c rerender_scalars" \
		-variable $this-scalars_transparency

	global $this-scalars_color_type
	make_labeled_radio $scalar.color \
	    "Scalar Coloring" "$this-c rerender_scalars" top 3 \
	    $this-scalars_color_type \
	    { {Default 0} {"Colormap Lookup" 1} {"RGB Conversion" 2} }

	sci_checkbutton $scalar.normalize \
		-text "Normalize before scaling" \
		-command "$this-c rerender_scalars" \
		-variable $this-scalars_normalize

	sci_checkbutton $scalar.smalldot \
		-text "Render glyphs below threshold" \
		-command "$this-c rerender_scalars" \
		-variable $this-scalars_small_is_dot

	make_labeled_radio $scalar.radio \
	    "Scalar Display Type" "$this-c rerender_scalars" top 2 \
	    $this-scalars_display_type \
	    {{Points Points} {Spheres Spheres} \
		 {Boxes Boxes} {Axes Axes}}
	
	pack $scalar.show_scalars $scalar.transparency \
	    $scalar.color $scalar.normalize $scalar.smalldot $scalar.radio \
	    -side top -pady 2 -fill y -anchor w

	global $this-ss_slider
	expscale $scalar.slide -label "ScalarScale" \
		-orient horizontal \
		-variable $this-scalars_scale
	set $this-ss_slider $scalar.slide
	bind $scalar.slide.scale <ButtonRelease> \
	    "$this-c scalars_scale; set $this-use_default_size 0"

	global $this-st_slider
	expscale $scalar.threshold -label "ScalarThreshold" \
		-orient horizontal \
		-variable $this-scalars_threshold
	set $this-st_slider $scalar.threshold
	bind $scalar.threshold.scale <ButtonRelease> \
	    "$this-c scalars_threshold; set $this-use_default_size 0"

	sci_labeledframe $scalar.resolution \
	    -labelpos nw -labeltext "Glyph Resolution"
	pack $scalar.resolution -side top -fill x -expand 1

	set res [$scalar.resolution childsite]
	sci_scale $res.scale -orient horizontal \
	    -variable $this-scalars_resolution \
	    -from 3 -to 20 -showvalue true -resolution 1
	bind $res.scale <ButtonRelease> "$this-c scalars_resolution"
	pack $res.scale -side top -pady 2 -fill both -expand 1
    }


    # Vector Tab
    method add_vector_tab {dof} {

	set vector [$dof.tabs add -label "Vectors" \
		-command "$this set_active_tab \"Vectors\""]
	sci_checkbutton $vector.show \
		-text "Show Vectors" \
		-command "$this-c toggle_display_vectors" \
		-variable $this-vectors_on

	sci_checkbutton $vector.transparent \
		-text "Enable Transparency" \
		-command "$this-c rerender_vectors" \
		-variable $this-vectors_transparency

	global $this-vectors_color_type
	make_labeled_radio $vector.color \
	    "Vector Coloring" "$this-c rerender_vectors" top 3 \
	    $this-vectors_color_type \
	    { {Default 0} {"Colormap Lookup" 1} {"RGB Conversion" 2} }

	sci_checkbutton $vector.normalize \
		-text "Normalize before scaling" \
		-command "$this-c rerender_vectors" \
		-variable $this-vectors_normalize

	sci_checkbutton $vector.bidirectional \
		-text "Render bidirectionally" \
		-command "$this-c rerender_vectors" \
		-variable $this-vectors_bidirectional

	sci_checkbutton $vector.smalldot \
		-text "Render glyphs below threshold" \
		-command "$this-c rerender_vectors" \
		-variable $this-vectors_small_is_dot

	make_labeled_radio $vector.radio \
	    "Vector Display Type" "$this-c rerender_vectors" top 4 \
	    $this-vectors_display_type \
	    {{Lines Lines} {Needles Needles} \
		 {Comets Comets} \
		 {Arrows Arrows} {Cones Cones} \
		 {Disks Disks} {Rings Rings} \
		 {Springs Springs}}
	
	pack $vector.show $vector.transparent \
	    $vector.color $vector.normalize \
	    $vector.bidirectional \
      $vector.smalldot \
	    $vector.radio -side top -pady 2 -fill y -anchor w

	global $this-vs_slider
	expscale $vector.slide -label "VectorScale" \
		-orient horizontal \
		-variable $this-vectors_scale
	set $this-vs_slider $vector.slide
	bind $vector.slide.scale <ButtonRelease> \
	    "$this-c vectors_scale; set $this-use_default_size 0"

	global $this-vt_slider
	expscale $vector.threshold -label "VectorThreshold" \
		-orient horizontal \
		-variable $this-vectors_threshold
	set $this-st_slider $vector.threshold
	bind $vector.threshold.scale <ButtonRelease> \
	    "$this-c vectors_threshold; set $this-use_default_size 0"

	sci_labeledframe $vector.resolution \
	    -labelpos nw -labeltext "Glyph Resolution"
	pack $vector.resolution -side top -fill x -expand 1

	set res [$vector.resolution childsite]
	sci_scale $res.scale -orient horizontal \
	    -variable $this-vectors_resolution \
	    -from 3 -to 20 -showvalue true -resolution 1
	bind $res.scale <ButtonRelease> "$this-c vectors_resolution"
	pack $res.scale -side top -fill both -expand 1
    }


    # Tensor Tab
    method add_tensor_tab {dof} {

	set tensor [$dof.tabs add -label "Tensors" \
		-command "$this set_active_tab \"Tensors\""]

	sci_checkbutton $tensor.show \
		-text "Show Tensors" \
		-command "$this-c toggle_display_tensors" \
		-variable $this-tensors_on

	sci_checkbutton $tensor.transparent \
		-command "$this-c rerender_tensors" \
		-text "Enable Transparency" \
		-variable $this-tensors_transparency

	global $this-tensors_color_type
	make_labeled_radio $tensor.color \
	    "Tensor Coloring" "$this-c rerender_tensors" top 3 \
	    $this-tensors_color_type \
	    { {Default 0} {"Colormap Lookup" 1} {"RGB Conversion" 2} }

	sci_checkbutton $tensor.normalize \
		-text "Normalize before scaling" \
		-command "$this-c rerender_tensors" \
		-variable $this-tensors_normalize

	sci_checkbutton $tensor.smalldot \
		-text "Render glyphs below threshold" \
		-command "$this-c rerender_tensors" \
		-variable $this-tensors_small_is_dot

	make_labeled_radio $tensor.radio \
	    "Tensor Display Type" "$this-c rerender_tensors" top 2 \
	    $this-tensors_display_type \
	    {{Boxes Boxes}  \
		 {"Colored Boxes" "Colored Boxes"}\
		 {Ellipsoids Ellipsoids} \
		 {Superquadrics Superquadrics}}
	
	pack $tensor.show $tensor.transparent \
	    $tensor.color $tensor.normalize \
      $tensor.smalldot \
	    $tensor.radio -side top -pady 2 -fill y -anchor w
	
	sci_labeledframe $tensor.emphasis \
	    -labelpos nw -labeltext "Superquadric Emphasis"
	pack $tensor.emphasis -side top -fill x -expand 1

	set emphasis [$tensor.emphasis childsite]

	sci_scale $emphasis.scale -orient horizontal \
	    -variable $this-tensors_emphasis -showvalue false \
	    -from 0.0 -to 1.0 -resolution 0.02

	bind $emphasis.scale <ButtonRelease> "$this-c tensors_resolution"
	pack $emphasis.scale -side top -fill both -expand 1


	global $this-ts_slider
	expscale $tensor.slide -label "TensorScale" \
		-orient horizontal \
		-variable $this-tensors_scale
	set $this-ts_slider $tensor.slide
	bind $tensor.slide.scale <ButtonRelease> \
	    "$this-c tensors_scale; set $this-use_default_size 0"

	global $this-tt_slider
	expscale $tensor.threshold -label "TensorThreshold" \
		-orient horizontal \
		-variable $this-tensors_threshold
	set $this-tt_slider $tensor.threshold
	bind $tensor.threshold.scale <ButtonRelease> \
	    "$this-c tensors_threshold; set $this-use_default_size 0"

	sci_labeledframe $tensor.resolution \
	    -labelpos nw -labeltext "Glyph Resolution"
	pack $tensor.resolution -side top -fill x -expand 1

	set res [$tensor.resolution childsite]
	sci_scale $res.scale -orient horizontal \
	    -variable $this-tensors_resolution \
	    -from 3 -to 20 -showvalue true -resolution 1
	bind $res.scale <ButtonRelease> "$this-c tensors_resolution"
	pack $res.scale -side top -fill both -expand 1
    }

    # Secondary Tab
    method add_secondary_tab {dof} {

	set secondary [$dof.tabs add -label "Secondary" \
		-command "$this set_active_tab \"Secondary\""]

 	sci_checkbutton $secondary.show \
 		-text "Use Secondary" \
		-command "$this-c rerender_all" \
 		-variable $this-secondary_on

	global $this-secondary_color_type
	make_labeled_radio $secondary.color_type \
	    "Use data for color " "$this-c rerender_all" top 3 \
	    $this-secondary_color_type \
	    { {Off 0} {"Colormap Lookup" 1} {"RGB Conversion" 2} }

 	sci_checkbutton $secondary.alpha \
 		-text "Use data for alpha mapping" \
		-command "$this-c rerender_all" \
 		-variable $this-secondary_alpha
 	sci_checkbutton $secondary.value \
 		-text "Use data for secondary glyph value" \
		-command "$this-c rerender_all" \
 		-variable $this-secondary_value

	make_labeled_radio $secondary.radio \
	    "Spring Display Type" "$this-c rerender_vectors" top 2 \
	    $this-secondary_display_type \
	    {{"Major Radius" "Major Radius"}
	     {"Minor Radius" "Minor Radius"} {Pitch Pitch} }

	pack $secondary.show $secondary.color_type \
	    $secondary.alpha $secondary.value \
	    $secondary.radio \
	    -side top -pady 2 -fill y -anchor w
	
	global $this-second_slider
	expscale $secondary.slide -label "SecondaryScale" \
		-orient horizontal \
		-variable $this-secondary_scale
	set $this-second_slider $secondary.slide
	bind $secondary.slide.scale <ButtonRelease> \
	    "$this-c rerender_all; set $this-use_default_size 0"
    }

    # Tertiary Tab
    method add_tertiary_tab {dof} {

	set tertiary [$dof.tabs add -label "Tertiary" \
		-command "$this set_active_tab \"Tertiary\""]

 	sci_checkbutton $tertiary.show \
 		-text "Use Tertiary" \
		-command "$this-c rerender_all" \
 		-variable $this-tertiary_on

	global $this-tertiary_color_type
	make_labeled_radio $tertiary.color_type \
	    "Use data for color " "$this-c rerender_all" top 3 \
	    $this-tertiary_color_type \
	    { {Off 0} {"Colormap Lookup" 1} {"RGB Conversion" 2} }

 	sci_checkbutton $tertiary.alpha \
 		-text "Use data for alpha mapping" \
		-command "$this-c rerender_all" \
 		-variable $this-tertiary_alpha
 	sci_checkbutton $tertiary.value \
 		-text "Use data for tertiary glyph value" \
		-command "$this-c rerender_all" \
 		-variable $this-tertiary_value

	make_labeled_radio $tertiary.radio \
	    "Spring Display Type" "$this-c rerender_vectors" top 2 \
	    $this-tertiary_display_type \
	    {{"Major Radius" "Major Radius"}
	     {"Minor Radius" "Minor Radius"} {Pitch Pitch} }

	pack $tertiary.show $tertiary.color_type \
	    $tertiary.alpha $tertiary.value \
	    $tertiary.radio \
	    -side top -pady 2 -fill y -anchor w
	
	global $this-third_slider
	expscale $tertiary.slide -label "TertiaryScale" \
		-orient horizontal \
		-variable $this-tertiary_scale
	set $this-third_slider $tertiary.slide
	bind $tertiary.slide.scale <ButtonRelease> \
	    "$this-c rerender_all; set $this-use_default_size 0"
    }

    method scalar_tab_changed {name1 name2 op} {
	global $this-scalars_has_data

	set window .ui[modname]
	if {[winfo exists $window]} {
	    set dof [$window.options.disp.frame_title childsite]	
	    if {[set $name1] == 1} { 
		add_scalar_tab $dof
		$dof.tabs view [set $this-active_tab]
	    } else {
		$dof.tabs delete "Scalars"
	    }
	}
    }

    method vector_tab_changed {name1 name2 op} {
	global $this-vectors_has_data

	set window .ui[modname]
	if {[winfo exists $window]} {
	    set dof [$window.options.disp.frame_title childsite]	
	    if {[set $name1] == 1} { 
		add_vector_tab $dof
		$dof.tabs view [set $this-active_tab]
	    } else {
		$dof.tabs delete "Vectors"
	    }
	}
    }

    method tensor_tab_changed {name1 name2 op} {
	global $this-tensors_has_data

	set window .ui[modname]
	if {[winfo exists $window]} {
	    set dof [$window.options.disp.frame_title childsite]	
	    if {[set $name1] == 1} { 
		add_tensor_tab $dof
		$dof.tabs view [set $this-active_tab]
	    } else {
		$dof.tabs delete "Tensors"
	    }
	}
    }

    method secondary_tab_changed {name1 name2 op} {
	global $this-secondary_has_data

	set window .ui[modname]
	if {[winfo exists $window]} {
	    set dof [$window.options.disp.frame_title childsite]	
	    if {[set $name1] == 1} { 
		add_secondary_tab $dof
		$dof.tabs view [set $this-active_tab]
	    } else {
		$dof.tabs delete "Secondary"
	    }
	}
    }

    method tertiary_tab_changed {name1 name2 op} {
	global $this-tertiary_has_data

	set window .ui[modname]
	if {[winfo exists $window]} {
	    set dof [$window.options.disp.frame_title childsite]	
	    if {[set $name1] == 1} { 
		add_tertiary_tab $dof
		$dof.tabs view [set $this-active_tab]
	    } else {
		$dof.tabs delete "Tertiary"
	    }
	}
    }

    method ui {} {
	set window .ui[modname]
	if {[winfo exists $window]} {
	    return
	}
	sci_toplevel $window
	#wm minsize $window 380 548

	#frame for all options to live
	sci_frame $window.options
 
	# node frame holds ui related to vert display (left side)
	sci_frame $window.options.disp -borderwidth 2
	pack $window.options.disp -padx 2 -pady 2 -side left \
		-fill both -expand 1

	# Display Options
	sci_labeledframe $window.options.disp.frame_title \
		-labelpos nw -labeltext "Display Options"
	set dof [$window.options.disp.frame_title childsite]

	sci_tabnotebook  $dof.tabs -height 490 -width 325 \
	    -raiseselect true 
	#label $window.options.disp.frame_title -text "Display Options"

	global $this-scalars_has_data
	global $this-vectors_has_data
	global $this-tensors_has_data
	global $this-secondary_has_data
	global $this-tertiary_has_data

	if {[set $this-scalars_has_data] == 1} {
	    add_scalar_tab $dof
	}

	if {[set $this-vectors_has_data] == 1} {
	    add_vector_tab $dof
	}

	if {[set $this-tensors_has_data] == 1} {
	    add_tensor_tab $dof
	}

	if {[set $this-secondary_has_data] == 1} {
	    add_secondary_tab $dof
	}

	if {[set $this-tertiary_has_data] == 1} {
	    add_tertiary_tab $dof
	}

	add_text_tab $dof

	global $this-active_tab
	global $this-interactive_mode
	# view the active tab
	if [catch "$dof.tabs view [set $this-active_tab]"] {
	    catch "$dof.tabs view 0"
	}

	$dof.tabs configure -tabpos "n"

	pack $dof.tabs -side top -fill x -expand yes -padx 2 -pady 2

	#pack notebook frame
	pack $window.options.disp.frame_title -side top -expand yes -fill x
	
	#add bottom frame for execute and dismiss buttons
	sci_frame $window.control -relief groove -borderwidth 2 -width 500
	sci_frame $window.def
	sci_frame $window.def.vals
	sci_frame $window.def.col
	sci_frame $window.def.col.f
	sci_frame $window.def.col.le

	pack $window.def.col $window.def.vals -side left -padx 10
	sci_label $window.def.col.le.approxl -text "PWL Approx Div:"
	sci_entry $window.def.col.le.approx -textvar $this-approx-div -width 4

	bind $window.def.col.le.approx <KeyRelease> "$this-c approx"



	addColorSelection $window.def.col.f "Default Color" \
	    $this-def_color "default_color_change"

	sci_button $window.def.vals.calcdefs -text "Calculate Default Size" \
		-command "$this-c calcdefs; set $this-use_default_size 1"
	sci_checkbutton $window.def.vals.use_defaults \
		-text "Use Default Size" \
		-variable $this-use_default_size

	pack $window.def.col.f $window.def.col.le -side top -pady 2 -anchor w
	pack $window.def.col.le.approxl $window.def.col.le.approx -side left
	pack $window.def.vals.use_defaults $window.def.vals.calcdefs \
	    -side top -pady 2

	sci_frame $window.fname -borderwidth 2
	sci_label $window.fname.label -text "Field Name"
	sci_entry $window.fname.entry -textvar $this-field_name
	sci_checkbutton $window.fname.override \
		-text "Override" \
		-command "$this-c rerender_all" \
		-variable $this-field_name_override

	TooltipMultiWidget "$window.fname.entry $window.fname.label" \
	    "Enter (optional) Field Name here.  The name will be displayed\nin the Viewer Window's list of Objects."
	
	pack $window.fname.label $window.fname.entry -side left
	pack $window.fname.override -side left -padx 6
	pack $window.fname -anchor w -padx 6 -pady 6

	# execute policy
	make_labeled_radio $window.control.exc_policy \
		"Execute Policy" "$this-c execute_policy" top 1 \
		$this-interactive_mode \
		{{"Interactively update" Interactive} \
		{"Execute button only" OnExecute}}

	#pack $res.scale -side top -fill both -expand 1

	pack $window.options -padx 2 -pady 2 -side top -fill x -expand 1
	#pack $window.resolution -padx 2 -pady 2 -side top -fill x -expand 1
	pack $window.def $window.control \
	    -padx 2 -pady 2 -side top

	pack $window.control.exc_policy -side top -fill both

	sci_frame $window.control.excdis -borderwidth 2
	pack $window.control.excdis -padx 4 -pady 4 -side top -fill both

	makeSciButtonPanel $window $window $this
	moveToCursor $window

	pack $window.control -padx 4 -pady 4 -side top -fill both
    }
}
