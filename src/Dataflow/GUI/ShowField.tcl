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


itcl::class SCIRun_Visualization_ShowField {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
	set name ShowField
	set_defaults
    }

    method set_defaults {} {

#	trace variable $this-nodes_scaleNV w "$this new_nodes_scale"
#	trace variable $this-edges_scaleNV w "$this new_edges_scale"
#	trace variable $this-active_tab w "$this switch_to_active_tab"

	# no C side component for these variables
	global $this-ns_slider
	set $this-ns_slider "not.set.yet"
	global $this-es_slider
	set $this-es_slider "not.set.yet"
    }

    method new_nodes_scale {a1 a2 a3} {
      global $this-ns_slider
      set val [set $this-nodes_scaleNV]
      if {$val != -0.0} {
        upvar $this-ns_slider ns_slider
        if {[info exists $this-ns_slider] && [winfo exists $ns_slider]} {
            $ns_slider newvalue $val
        } else {
            set $this-nodes_scale [set $this-nodes_scaleNV]
        }
        set $this-nodes_scaleNV -0.0
      }
    }
    
    
    method new_edges_scale {a1 a2 a3} {
      global $this-es_slider
      set val [set $this-edges_scaleNV]
      if {$val != -0.0} {
        upvar $this-es_slider es_slider
        if {[info exists $this-es_slider] && [winfo exists $es_slider]} {
            $es_slider newvalue $val
        } else {
            set $this-edges_scale [set $this-edges_scaleNV]
        }
        set $this-edges_scaleNV -0.0
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
	set window .ui[modname]
	if {[winfo exists $window]} {
	    set dof [$window.options.disp.frame_title childsite]
	    $dof.tabs view [set $this-active_tab]
	}
    }

    # Nodes Tab
    method add_nodes_tab {dof inserting} {
	
	if {$inserting} {
	    set node [$dof.tabs insert 0 -label "Nodes" \
			  -command "$this set_active_tab \"Nodes\""]
	} else {
	    set node [$dof.tabs add -label "Nodes" \
			  -command "$this set_active_tab \"Nodes\""]
	}
	
	sci_checkbutton $node.show \
		-text "Show Nodes" \
		-command "$this-c toggle_display_nodes" \
		-variable $this-nodes_on
	sci_checkbutton $node.transparency \
		-text "Enable Transparency" \
		-command "$this-c rerender_nodes" \
		-variable $this-nodes_transparency

	global $this-nodes_color_type
	
	make_labeled_radio $node.color \
	    "Node Coloring " "$this-c rerender_nodes" top 1\
	    $this-nodes_color_type \
	    { {Default 0} {"Colormap Lookup" 1} \
		  {"Conversion to RGB" 2} }

	global $this-nodes_display_type
	
	make_labeled_radio $node.radio \
	    "Node Display Type" "$this-c rerender_nodes" top 1 \
	    $this-nodes_display_type \
	    { {Points Points} {Spheres Spheres} }

	pack $node.show $node.transparency $node.color $node.radio \
	    -fill y -anchor w

	global $this-ns_slider
	expscale $node.slide -label SphereScale \
	    -orient horizontal \
	    -variable $this-nodes_scale

	set $this-ns_slider $node.slide

	bind $node.slide.scale <ButtonRelease> \
	    "$this-c nodes_scale; set $this-use_default_size 0"

	sci_labeledframe $node.resolution \
	    -labelpos nw -labeltext "Resolution"
	pack $node.resolution -side top -fill x -expand 1

	set res [$node.resolution childsite]
	sci_scale $res.scale -orient horizontal -variable $this-nodes_resolution \
	    -from 3 -to 20 -showvalue true -resolution 1
	bind $res.scale <ButtonRelease> "$this-c rerender_nodes"
	pack $res.scale -side top -fill both -expand 1
    }

    # Edges Tab
    method add_edges_tab {dof} {

	set edge [$dof.tabs add -label "Edges" \
		-command "$this set_active_tab \"Edges\""]

	sci_checkbutton $edge.show \
		-text "Show Edges" \
		-command "$this-c toggle_display_edges" \
		-variable $this-edges_on
	sci_checkbutton $edge.transparency \
		-text "Enable Transparency" \
		-command "$this-c rerender_edges" \
		-variable $this-edges_transparency

	global $this-edges_color_type
	
	make_labeled_radio $edge.color \
	    "Edge Coloring " "$this-c rerender_edges" top 1\
	    $this-edges_color_type \
	    { {Default 0} {"Colormap Lookup" 1} \
		  {"Conversion to RGB" 2} }

	make_labeled_radio $edge.radio \
		"Edge Display Type" "$this-c rerender_edges" top 1\
		$this-edges_display_type {{Lines Lines} {Cylinders Cylinders}}

	pack $edge.show $edge.transparency $edge.color $edge.radio \
		-side top -fill y -anchor w
	global $this-es_slider
	expscale $edge.slide -label CylinderScale \
		-orient horizontal \
		-variable $this-edges_scale
	set $this-es_slider $edge.slide

	bind $edge.slide.scale <ButtonRelease> \
	    "$this-c edges_scale; set $this-use_default_size 0"

	sci_labeledframe $edge.resolution \
	    -labelpos nw -labeltext "Resolution"
	pack $edge.resolution -side top -fill x -expand 1

	set res [$edge.resolution childsite]
	sci_scale $res.scale -orient horizontal -variable $this-edges_resolution \
	    -from 3 -to 20 -showvalue true -resolution 1
	bind $res.scale <ButtonRelease> "$this-c rerender_edges"
	pack $res.scale -side top -fill both -expand 1
    }

    # Faces Tab
    method add_faces_tab {dof} {
	set face [$dof.tabs add -label "Faces" \
		-command "$this set_active_tab \"Faces\""]
	sci_checkbutton $face.show \
		-text "Show Faces" \
		-command "$this-c toggle_display_faces" \
		-variable $this-faces_on
	sci_checkbutton $face.transparency \
		-text "Enable Transparency" \
		-command "$this-c rerender_faces" \
		-variable $this-faces_transparency

	global $this-faces_color_type
	make_labeled_radio $face.color \
	    "Face Coloring " "$this-c rerender_faces" top 1\
	    $this-faces_color_type \
	    { {Default 0} {"Colormap Lookup" 1} \
		  {"Conversion to RGB" 2} }

	sci_label $face.blank -text "  "

	sci_checkbutton $face.normals \
		-text "Use Face Normals" \
		-command "$this-c rerender_faces" \
		-variable $this-faces_normals
 	sci_checkbutton $face.texture \
 	        -text "Render Images as a texture (Colormap Only)" \
 	        -command "$this-c rerender_faces" \
 	        -variable $this-faces_usetexture

	pack $face.show $face.transparency $face.color \
	    $face.blank $face.normals $face.texture \
	    -side top -fill y -anchor w
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
	    "Text Size:" "$this-c rerender_text" left 5\
	    $this-text_fontsize \
	    {{"XS" 0} {"S" 1} {"M" 2} {"L" 3} {"XL" 4}}

	pack $text.show.data $text.show.nodes $text.show.edges \
	    $text.show.faces $text.show.cells \
	    -side top -fill y -anchor w
	
	global $this-text_color_type
	make_labeled_radio $text.color \
	    "Text Coloring " "$this-c rerender_text" top 1\
	    $this-text_color_type \
	    { {Text 0} {"Colormap Lookup" 1} \
		  {"Conversion to RGB" 2} }

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
	    $text.precision -side top -fill y -anchor w
    }

    method ui {} {
	set window .ui[modname]
	if {[winfo exists $window]} {
	    return
	}
	sci_toplevel $window
	#wm minsize $window 380 548

	#sci_frame for all options to live
	sci_frame $window.options
 
	# node sci_frame holds ui related to vert display (left side)
	sci_frame $window.options.disp -borderwidth 2
	pack $window.options.disp -padx 2 -pady 2 -side left \
		-fill both -expand 1

	# Display Options
	sci_labeledframe $window.options.disp.frame_title \
		-labelpos nw -labeltext "Display Options"
	set dof [$window.options.disp.frame_title childsite]

	sci_tabnotebook  $dof.tabs -height 420 -width 275 \
	    -raiseselect true 
	#label $window.options.disp.frame_title -text "Display Options"

	add_nodes_tab $dof 0
	add_edges_tab $dof
	add_faces_tab $dof
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
	
	#add bottom sci_frame for execute and dismiss buttons
	sci_frame $window.control -relief groove -borderwidth 2 -width 500
	sci_frame $window.def
	sci_frame $window.def.vals
	sci_frame $window.def.col
	sci_frame $window.def.col.f
	sci_frame $window.def.col.le

	pack $window.def.col $window.def.vals -side left -padx 10
	sci_label $window.def.col.le.approxl -text "PWL Approx Div:"
	sci_entry $window.def.col.le.approx -textvar $this-approx_div -width 4

	bind $window.def.col.le.approx <KeyRelease> "$this-c approx"

	addColorSelection $window.def.col.f "Default Color" \
	    $this-def_color "default_color_change"

	sci_checkbutton $window.def.vals.use_defaults \
		-text "Use Default Size" \
		-variable $this-use_default_size
	sci_button $window.def.vals.calcdefs -text "Calculate Default Size" \
		-command "$this-c calcdefs; set $this-use_default_size 1"

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


	## Cylinder and Sphere Resolution
	#iwidgets::labeledframe $window.resolution \
	#	-labelpos nw -labeltext "Cylinder and Sphere Resolution"
	#set res [$window.resolution childsite]
	#
	#scale $res.scale -orient horizontal -variable $this-resolution \
	#	-from 3 -to 20 -showvalue true -resolution 1
	#
	#bind $res.scale <ButtonRelease> "$this-c resolution_scale"

	# execute policy
	make_labeled_radio $window.control.exc_policy \
		"Execute Policy" "$this-c execute_policy" top 1\
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

	pack $window.control -padx 4 -pady 4 -side top -fill both

	makeSciButtonPanel $window $window $this
	moveToCursor $window

    }
}
