#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2011 Scientific Computing and Imaging Institute,
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

itcl::class BioPSE_Forward_BuildBEMatrix {
  inherit Module

  constructor { {args ""} } {
    eval configure $args
    set name BuildBEMatrix

    set_defaults
  }

  method set_defaults {} {
    initGlobal $this-inputs-listbox ""

    initGlobal $this-input_fields []
    initGlobal $this-input_field_types []
    initGlobal $this-input_field_inside_conductivity []
    initGlobal $this-input_field_outside_conductivity []
    initGlobal $this-input_field_surface_type []

    initGlobal $this-selected_field_index 0
    initGlobal $this-selected_inside_conductivity ""
    initGlobal $this-selected_outside_conductivity ""
    initGlobal $this-selected_surface_type 0

    initGlobal $this-null-frame "null"

    set $this-top_frame 0
  }

  method clear_lists {} {
    set $this-input_fields []
    set $this-input_field_types []
    set $this-input_field_inside_conductivity []
    set $this-input_field_outside_conductivity []
    set $this-input_field_surface_type []
  }

  method update_fields {} {
    global $this-input_field_list
    global $this-field_type_list
    global $this-inside_cond_list
    global $this-outside_cond_list
    global $this-surface_type_list
    global $this-input_fields
    global $this-top_frame
    global $this-selected_field
    global $this-selected_field_type
    global $this-selected_surface_type

    clear_lists

    set field_list [split [set $this-input_field_list]]
    set $this-input_fields $field_list

    set field_type_list [split [set $this-field_type_list]]
    set $this-input_field_types $field_type_list

    set inside_cond_list [split [set $this-inside_cond_list]]
    set $this-input_field_inside_conductivity $inside_cond_list

    set outside_cond_list [split [set $this-outside_cond_list]]
    set $this-input_field_outside_conductivity $outside_cond_list

    set surface_type_list [split [set $this-surface_type_list]]
    set $this-input_field_surface_type $surface_type_list

    if {[set $this-top_frame] == 0} {
      return;
    }

    set top_frame [set $this-top_frame]
    set index 0

    if {[ expr [ llength [set $this-input_fields] ] ] > 0} {
      if {[winfo exists $top_frame]} {
        foreach f $field_list {
	  make_frame $top_frame.$f [lindex [set $this-input_field_types] $index]
          incr index
        }

        # TODO: if a selected field is set, but other variables are not, then this part will break
        # Is there a better Tk solution?
        if [expr [string equal [set $this-selected_field] ""] == 1] {
          set $this-selected_field_index 0
          set $this-selected_field [lindex [set $this-input_fields] 0]
          set $this-selected_field_type [lindex [set $this-input_field_types] 0]
          set $this-selected_inside_conductivity [lindex [set $this-input_field_inside_conductivity] 0]
          set $this-selected_outside_conductivity [lindex [set $this-input_field_outside_conductivity] 0]
          set $this-selected_surface_type [lindex [set $this-input_field_surface_type] 0]

          set field [set $this-selected_field]
          pack $top_frame.$field -anchor nw -side left -fill both -expand yes -padx 2 -pady 2
        }
      } else {
        # TODO: needed?
        set fname [set $this-null-frame]
        make_null_frame $top_frame.$fname
        pack $top_frame.$fname -anchor nw -side left -fill both -expand yes -padx 2 -pady 2
      }
    }
  }

  method ui {} {
    set w .ui[modname]

    if {[winfo exists $w]} {
      return
    }
    global $this-selected_field
    global $this-input_field_list

    sci_toplevel $w
    wm minsize $w 300 100

    # TODO: maybe grid would be better than pack?
    sci_frame $w.f
    set $this-top_frame $w.f
    pack $w.f \
      -anchor nw -side top -fill both -expand yes -padx 2 -pady 2

    sci_frame $w.f.inputs
    pack $w.f.inputs \
      -anchor nw -side left -fill both -expand yes -padx 2 -pady 2

    sci_scrolledlistbox $w.f.inputs.listbox -selectioncommand [format "%s choose_input" $this]
    set $this-inputs-listbox $w.f.inputs.listbox
    $w.f.inputs.listbox component listbox configure -listvariable $this-input_field_list -selectmode browse
    pack $w.f.inputs.listbox \
      -anchor nw -side top -fill both -expand yes -padx 2 -pady 2

    update_fields

    makeSciButtonPanel $w $w $this
    moveToCursor $w
  }

  method choose_input {} {
    global $this-selected_field
    global $this-selected_field_type
    global $this-selected_surface_type
    global $this-input_field_list
    global $this-inputs-listbox
    global $this-input-type
    global $this-top_frame

    set top_frame [set $this-top_frame]
    set field [set $this-selected_field]

    if {[winfo exists $top_frame.$field] && [winfo ismapped $top_frame.$field]} {
      pack forget $top_frame.$field
    }

    # update from previously visible frame
    update_properties

    set inputnum [[set $this-inputs-listbox] curselection]
    set $this-selected_field_index $inputnum

    if [expr [string equal $inputnum ""] == 0] {
      set $this-selected_field  [lindex [set $this-input_field_list] $inputnum]
      set $this-selected_field_type [lindex [set $this-input_field_types] $inputnum]
      set $this-selected_inside_conductivity [lindex [set $this-input_field_inside_conductivity] $inputnum]
      set $this-selected_outside_conductivity [lindex [set $this-input_field_outside_conductivity] $inputnum]
      set $this-selected_surface_type [lindex [set $this-input_field_surface_type] $inputnum]

      set field [set $this-selected_field]
      pack $top_frame.$field -anchor nw -side left -fill both -expand yes -padx 2 -pady 2
    }
  }

  method make_null_frame {frame} {
    if {[winfo exists $frame]} {
      if {[winfo ismapped $frame]} {
        pack forget $frame
      }
      destroy $frame
    }

    sci_frame $frame
  }

  method make_frame {frame field_type} {
    if {[winfo exists $frame]} {
      if {[winfo ismapped $frame]} {
        pack forget $frame
      }
      destroy $frame
    }

    global $this-selected_field
    global $this-selected_field_type

    sci_frame $frame

    sci_frame $frame.name_field
    pack $frame.name_field \
      -anchor nw -side left -expand yes -padx 1p -pady 1p
 
    sci_frame $frame.input_type_field
    pack $frame.input_type_field \
      -anchor nw -side left -expand yes -padx 1p -pady 1p

    sci_frame $frame.module_properties
    pack $frame.module_properties \
      -anchor nw -side left -expand yes -padx 1p -pady 1p

    sci_label $frame.name_field.namelabel -text "Name:"
    pack $frame.name_field.namelabel \
      -anchor nw -side left -padx 1p -pady 1p

    sci_entry $frame.name_field.name \
      -textvariable $this-selected_field -state readonly
    pack $frame.name_field.name \
      -anchor nw -side right -padx 1p -pady 1p

    # could use sci_labeledframes instead
    sci_label $frame.input_type_field.label -text "Field Type:"
    pack $frame.input_type_field.label \
      -anchor nw -side top -padx 2p -pady 2p
 
    sci_entry $frame.input_type_field.field_type \
      -textvariable $this-selected_field_type -state readonly
    pack $frame.input_type_field.field_type \
      -anchor nw -side top -padx 2p -pady 2p

    sci_frame $frame.input_type_field.type_subframe
    pack $frame.input_type_field.type_subframe \
      -anchor nw -side top -expand yes -padx 1p -pady 1p

    sci_radiobutton $frame.input_type_field.type_subframe.measurement \
      -variable $this-selected_surface_type -text "Measurement (Neumann Boundary Condition)" -value 1 \
      -command [format "%s update_properties" $this]
    pack $frame.input_type_field.type_subframe.measurement \
      -anchor nw -side left -padx 2p -pady 2p
    sci_radiobutton $frame.input_type_field.type_subframe.source \
      -variable $this-selected_surface_type -text "Source (Dirichlet Boundary Condition)" -value 0 \
      -command [format "%s update_properties" $this]
    pack $frame.input_type_field.type_subframe.source \
      -anchor nw -side right -padx 2p -pady 2p

    make_entry $frame.module_properties.inside "Inside Conductivity:" $this-selected_inside_conductivity [format "%s update_properties" $this]
    make_entry $frame.module_properties.outside "Outside Conductivity:" $this-selected_outside_conductivity [format "%s update_properties" $this]
    pack $frame.module_properties.inside $frame.module_properties.outside \
      -anchor nw -side top -padx 2p -pady 2p

    sci_frame $frame.update_fields

    sci_button $frame.update_fields.update \
      -text "Update Field Properties" -command [format "%s update_properties" $this]
    pack $frame.update_fields.update \
      -anchor nw -side top -padx 2p -pady 2p

    pack $frame.name_field $frame.input_type_field $frame.input_type_field \
      $frame.module_properties $frame.update_fields \
      -anchor nw -side top -expand yes -padx 1p -pady 1p
  }

  method make_entry {w text textvar command} {
    if {[winfo exists $w]} {
      return
    }
    sci_frame $w
    sci_label $w.l -text "$text"
    pack $w.l -anchor nw -side left -padx 2p -pady 2p
    global $textvar
    sci_entry $w.e -textvariable $textvar
    bind $w.e <Return> $command
    pack $w.e -anchor nw -side right -padx 2p -pady 2p
  }

  method update_properties {} {
    # TODO: is there any way to validate input in Tk? (i.e. check for numerical input)
    global $this-selected_field_index

    set index [set $this-selected_field_index]

    # inside conductivity property
    lset $this-input_field_inside_conductivity $index [set $this-selected_inside_conductivity]
    set $this-inside_cond_list [join [set $this-input_field_inside_conductivity] ]

    # outside conductivity property
    lset $this-input_field_outside_conductivity $index [set $this-selected_outside_conductivity]
    set $this-outside_cond_list [join [set $this-input_field_outside_conductivity] ]

    # surface type property (measurement or source)
    lset $this-input_field_surface_type $index [set $this-selected_surface_type]
    set $this-surface_type_list [join [set $this-input_field_surface_type] ]
  }

  # TODO: may not be needed
  method activate { w } {
    $w configure -state normal -foreground black
  }

  # TODO: may not be needed
  method deactivate { w } {
    $w configure -state disabled -foreground darkgrey
  }
}