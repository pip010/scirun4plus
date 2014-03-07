#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2012 Scientific Computing and Imaging Institute,
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


catch {rename BioPSE_Inverse_SolveInverseProblemWithTikhonov ""}

itcl::class BioPSE_Inverse_SolveInverseProblemWithTikhonov {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name SolveInverseProblemWithTikhonov
        set_defaults
    }

    method set_defaults {} {
      global $this-lambda_lc
      global $this-lambda_index

      # lambda_lc is the lambda corner value set when the lcurve is plotted
      # using plot_graph, does not have a C++ counterpart
      set $this-lambda_lc 0.0
      set $this-lambda_index 0
      set $this-graph_needs_updating 0
      set $this-epsilon 1e-9

      set $this-lcurve []
      set $this-lcorner []
    }

    method entervalue {} {
      set w .ui[modname]
      set color "#505050"
      $w.ent.f2.l1 configure -foreground black
      $w.ent.f2.e1 configure -state normal -foreground black
      $w.sld.sld.f3.s configure -state disabled -foreground $color
      $w.lc.lam.l configure -foreground $color
      $w.lc.lam.e configure -state disabled -foreground $color
      clear_graph
      $w.sld.sld.f2.lvalue configure -state disabled -foreground $color
      $w.sld.sld.f2.l1 configure -foreground $color
      $w.lc.f2.e1 configure -state disabled -foreground $color
      $w.lc.f2.l1 configure -foreground $color
    }

    method useslider {} {
      set w .ui[modname]
      set color "#505050"
      $this calc_log [set $this-lambda_from_scale]
      $w.sld.sld.f3.s configure -state normal -foreground black
      $w.ent.f2.e1 configure -state disabled -foreground $color
      $w.ent.f2.l1 configure -foreground $color
      $w.lc.lam.l configure -foreground $color
      $w.lc.lam.e configure -state disabled -foreground $color
      clear_graph
      $w.sld.sld.f2.lvalue configure -state disabled -foreground black
      $w.sld.sld.f2.l1 configure -foreground black
      $w.lc.f2.e1 configure -state disabled -foreground $color
      $w.lc.f2.l1 configure -foreground $color
    }

    method uselcurve {} {
      set w .ui[modname]
      set color "#505050"
      $w.sld.sld.f3.s configure -state disabled -foreground $color
      $w.ent.f2.l1 configure -foreground $color
      $w.ent.f2.e1 configure -state disabled -foreground $color
      $w.lc.lam.l configure -foreground black
      $w.lc.lam.e configure -state readonly -foreground black
      clear_graph
      $w.sld.sld.f2.lvalue configure -state disabled -foreground $color
      $w.sld.sld.f2.l1 configure -foreground $color
      $w.lc.f2.e1 configure -state normal -foreground black
      $w.lc.f2.l1 configure -foreground black

      set $this-graph_needs_updating 1
      $this-c updategraph [set $this-lambda_lc] [set $this-lambda_index]
    }

    method execrunmode {} {
      set w .ui[modname]
      $this-c needexecute
    }

    ##############
    method plot_graph {lcurv lcor lam lam_index} {
        set w .ui[modname]
        if {![winfo exists $w]} {
          return
        }

        global $this-lambda_lc
        global $this-lambda_index
        global $this-lcurve
        global $this-lcorner

        # necessary, otherwise GUI crashes with an 'Illegal instruction' error
        if { ( ! [set $this-graph_needs_updating] ) &&
             ( ! [ expr {abs($lam - [set $this-lambda_min])} ] < [set $this-epsilon] ) &&
             ( [ expr {abs($lam - [set $this-lambda_lc])} ] < [set $this-epsilon] || $lam_index == [set $this-lambda_index] ) } {
          set $this-lambda_index $lam_index
          return
        }

        set $this-graph_needs_updating 0
        set $this-lambda_lc $lam
        set $this-lambda_index $lam_index
        set $this-lcurve $lcurv
        set $this-lcorner $lcor

        $w.lc.data.g element delete LCurve
        $w.lc.data.g element delete LCorner
        $w.lc.data.g element create LCurve -data $lcurv -color blue -symbol ""
        $w.lc.data.g element create LCorner -data $lcor -color green -symbol ""
    }

    ##############
    method clear_graph {} {
        set w .ui[modname]
        if {![winfo exists $w]} {
            return
        }

        set color "#505050"
        $w.lc.data.g element delete LCurve
        $w.lc.data.g element delete LCorner
        $w.lc.data.g element create LCurve -data "1 0"
        $w.lc.data.g element create LCorner -data "1 0"
    }

    ##############
    method calc_log {slider_value} {
        global $this-log_val
        global $this-lambda_min
        global $this-lambda_max

        set $this-log_val [expr {[set $this-lambda_min]*pow(10,$slider_value/5*log10([set $this-lambda_max]/[set $this-lambda_min]))}]
    }

    ##############
    
    # TODO: replace slider with mouse interaction - interaction with graph is crude and not very responsive
    # method update_lambda_corner {slider_value} {
    #     set w .ui[modname]
    #     if {![winfo exists $w]} {
    #         return
    #     }
    #     global $this-lambda_min
    #     global $this-lambda_max
    #     global $this-lambda_lc
    #     global $this-lambda_index
    #     global $this-lcurve
    #     global $this-lcorner
    #     if { [ expr {abs($slider_value - [set $this-lambda_lc])} ] < [set $this-epsilon] } {
    #       return
    #     }
    #     if {$slider_value < [set $this-lambda_min]} {
    #       set slider_value [set $this-lambda_min]
    #     }
    #     set $this-lambda_lc $slider_value
    #     if {[expr [ llength [set $this-lcurve] ] ] > 0 && [expr [ llength [set $this-lcorner] ] ] > 0} {
    #       $this-c updategraph $slider_value [set $this-lambda_index]
    #     }
    # }

    # method update_lambda_entry {} {
    #   set $this-lambda_from_textentry [set $this-lambda_lc]
    #   set $this-reg_method "single"
    #   eval "$this entervalue"
    # }

    ##############	
    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return;
        }

        global $this-have_ui
        global medium_font_mod

        set $this-have_ui 1
        set color "#505050"

        ###################################
        #THIS PART IS FOR MAKING THE TITLE
        ###################################
        sci_toplevel $w
        wm minsize $w 750 20

        sci_frame $w.warning_frame -relief ridge -border 2
        pack $w.warning_frame -side top -padx 2 -pady 2 -expand yes -fill both
        sci_label $w.warning_frame.label -font $medium_font_mod \
          -text "NOTE: See documentation (press ? below) for explanation of options."
        pack $w.warning_frame.label -side top -anchor w -padx 2 -pady 4 -fill x
  
        # TODO: use a grid layout in Qt

        sci_frame $w.f7 -relief ridge -border 2
        pack $w.f7 -padx 2 -pady 2 -expand yes -fill both

        sci_label $w.f7.warning_label -relief sunken -border 1 \
          -text "The user is responsible for ensuring that constraint matrices sent as inputs to the module match the choices below."
        pack $w.f7.warning_label -side top -anchor w -padx 2 -pady 4 -fill x

        sci_label $w.f7.label -text "Regularization formulation selection:"
        pack $w.f7.label -side top -anchor w -padx 2

        sci_frame $w.f7.cases
        pack $w.f7.cases -padx 20 -pady 2 -expand yes -fill x

        sci_radiobutton $w.f7.cases.automatic -text "automatic" \
          -variable $this-tik_cases -value 0
        pack $w.f7.cases.automatic -side top -anchor w -padx 2

        sci_radiobutton $w.f7.cases.underdetermined -text "underdetermined " \
          -variable $this-tik_cases -value 1
        pack $w.f7.cases.underdetermined -side top -anchor w -padx 2

        sci_radiobutton $w.f7.cases.overdetermined  -text "overdetermined " \
          -variable $this-tik_cases -value 2
        pack $w.f7.cases.overdetermined -side top -anchor w -padx 2 

        sci_label $w.f7.label2 -text "Description of constraint matrices:"
        pack $w.f7.label2 -side top -anchor w -padx 2

        sci_frame $w.f7.subcases
        pack $w.f7.subcases -padx 20 -pady 2 -expand yes -fill x

        sci_label $w.f7.subcases.label1 -text "Solution constraint"
        pack $w.f7.subcases.label1 -side top -anchor w -padx 2

        sci_frame $w.f7.subcases.s1
        pack $w.f7.subcases.s1 -padx 20 -pady 2 -expand yes -fill x

        sci_radiobutton $w.f7.subcases.s1.solution_constraint -state normal \
          -text "Solution constraint matrix" -variable $this-tik_solution_subcases -value 0
        sci_radiobutton $w.f7.subcases.s1.solution_constraint_squared -state normal \
          -text "Squared solution constraint matrix" -variable $this-tik_solution_subcases -value 1
        pack $w.f7.subcases.s1.solution_constraint $w.f7.subcases.s1.solution_constraint_squared -side left -anchor w -padx 2

        sci_label $w.f7.subcases.label2 -text "Residual constraint"
        pack $w.f7.subcases.label2 -side top -anchor w -padx 2

        sci_frame $w.f7.subcases.s2
        pack $w.f7.subcases.s2 -padx 20 -pady 2 -expand yes -fill x

        sci_radiobutton $w.f7.subcases.s2.residual_constraint -state normal \
          -text "Residual constraint matrix" -variable $this-tik_residual_subcases -value 0
        sci_radiobutton $w.f7.subcases.s2.residual_constraint_squared -state normal \
          -text "Squared residual constraint matrix" -variable $this-tik_residual_subcases -value 1
        pack $w.f7.subcases.s2.residual_constraint $w.f7.subcases.s2.residual_constraint_squared -side left -anchor w -padx 2

        sci_frame $w.titlefr -relief ridge -border 2
        sci_label $w.titlefr.l -text "Select Method for Lambda"
        pack $w.titlefr -side top -padx 2 -pady 2 -expand yes -fill x
        pack $w.titlefr.l -side top -fill x

        global $this-lambda_from_textentry
        global $this-lambda_from_scale
        global $this-lambda_min
        global $this-lambda_max
        global $this-lambda_resolution
        global $this-reg_method

        #######################
        # Entry radio-button
        #######################
        sci_frame $w.ent -relief ridge -border 2
        pack $w.ent -side top -padx 2 -pady 2 -expand yes -fill both

        sci_frame $w.ent.f1 -relief flat
        pack $w.ent.f1 -side top -padx 2 -pady 4 -expand yes -fill both

        sci_radiobutton $w.ent.f1.b -text "Enter value" \
          -variable "$this-reg_method" -value single \
          -command "$this entervalue"

        pack $w.ent.f1.b -side left

        sci_frame $w.ent.f2 -relief flat
        pack $w.ent.f2 -side top -padx 2 -pady 2 -expand yes -fill both
        sci_label $w.ent.f2.l1 -text "Lambda: "
        sci_entry $w.ent.f2.e1 -textvariable $this-lambda_from_textentry
        pack $w.ent.f2.l1 $w.ent.f2.e1 -side left -expand yes \
          -fill x -padx 2 -pady 2
        #bind $w.ent.f2.e1 <Return> $n

        #########################
        # Slider radio-button
        #########################
        sci_frame $w.sld -relief flat
        pack $w.sld -side top -expand yes -fill x -padx 2 -pady 2
	
        sci_frame $w.sld.sld -relief ridge -border 2
        pack $w.sld.sld -side left -expand yes -fill x -padx 2 -pady 2

        sci_frame $w.sld.rng -relief ridge -border 2
        pack $w.sld.rng -side right  -fill y -padx 2 -pady 2
        
        sci_frame $w.sld.rng.f1 -relief flat
        pack  $w.sld.rng.f1 -side top -fill x
    
        # TODO: put a refresh button for From/Step Size/To text entry boxes...
        # or for Qt migration, refresh on focus out

        sci_label $w.sld.rng.f1.l1 -text "Lambda Range"	
        pack  $w.sld.rng.f1.l1 -side top -fill x

        sci_frame $w.sld.rng.f2 -relief flat
        pack  $w.sld.rng.f2 -side top -fill x -padx 2 -pady 2

        sci_entry $w.sld.rng.f2.e1 -textvariable $this-lambda_min
        sci_label $w.sld.rng.f2.l1 -text "From "
        pack $w.sld.rng.f2.e1  $w.sld.rng.f2.l1 -side right

        # inserted out of order...
        sci_frame $w.sld.rng.f4 -relief flat
        pack  $w.sld.rng.f4 -side top -fill x -padx 2 -pady 2

        sci_entry $w.sld.rng.f4.e1 -textvariable $this-lambda_resolution
        sci_label $w.sld.rng.f4.l1 -text "Step Size "
        pack $w.sld.rng.f4.e1  $w.sld.rng.f4.l1 -side right
        
        sci_frame $w.sld.rng.f3 -relief flat
        pack  $w.sld.rng.f3 -side top -fill x -padx 2 -pady 2

        sci_entry $w.sld.rng.f3.e1 -textvariable $this-lambda_max
        sci_label $w.sld.rng.f3.l1 -text "  To  "
        pack $w.sld.rng.f3.e1  $w.sld.rng.f3.l1 -side right

        sci_frame $w.sld.sld.f1 -relief flat
        pack $w.sld.sld.f1 -side top -expand yes -fill x -padx 2 -pady 4

        sci_radiobutton $w.sld.sld.f1.b -text "Choose using slider" \
          -variable "$this-reg_method" -value slider \
          -command "$this useslider"
        pack $w.sld.sld.f1.b -side left

        sci_frame $w.sld.sld.f2 -relief flat
        pack $w.sld.sld.f2 -side top -expand yes -fill x -padx 2 -pady 2
        
        sci_entry $w.sld.sld.f2.lvalue -textvariable "$this-log_val" -state disabled       
        sci_label $w.sld.sld.f2.l1 -text "Lambda: "
        pack $w.sld.sld.f2.l1 $w.sld.sld.f2.lvalue -side left
         
        sci_frame $w.sld.sld.f3 -relief flat
              pack $w.sld.sld.f3 -side top -expand yes -fill x -padx 2 -pady 2

        sci_scale $w.sld.sld.f3.s -from 0.0 -to 5.0 \
          -resolution 0.01 -orient horizontal \
          -variable $this-lambda_from_scale -showvalue false \
          -command "$this calc_log" 
        pack $w.sld.sld.f3.s -side left -expand yes -fill x -padx 2 -pady 2

        #######################
        # LCurve radio-button
        #######################

        sci_frame $w.lc -relief ridge -border 2
        pack $w.lc -side top -padx 2 -pady 2 -expand yes -fill both
	
        sci_frame $w.lc.f1 -relief flat
        pack $w.lc.f1 -side top -padx 2 -pady 2 -expand yes -fill x

        sci_radiobutton $w.lc.f1.b -text "L-curve             " \
          -variable "$this-reg_method" -value lcurve \
          -command "$this uselcurve"
        pack $w.lc.f1.b -side left
	
        sci_frame $w.lc.f2 -relief flat
        pack $w.lc.f2 -side top -padx 2 -pady 4 -expand yes -fill both

        sci_label $w.lc.f2.l1 -text "Number of Points: "
        sci_entry $w.lc.f2.e1 -textvariable $this-lambda_num 
        pack $w.lc.f2.e1 $w.lc.f2.l1 -side right	

        sci_frame $w.lc.data -relief groove -borderwidth 2
        blt::graph $w.lc.data.g -height 250 \
                -plotbackground #CCCCFF
        $w.lc.data.g element create LCurve -data "1 0" -color blue -symbol ""
        $w.lc.data.g element create LCorner -data "1 0" -color green -symbol ""
        $w.lc.data.g yaxis configure -title "log || Rx ||"
        $w.lc.data.g xaxis configure -title "log || Ax - y ||"

        pack $w.lc.data -side top -fill x -padx 2 -pady 2
        pack $w.lc.data.g -side top -fill x

        # TODO: replace slider with mouse interaction - interaction with graph is crude and not very responsive.
        # Mouse pointer should be able to pick point, then point should be interpolated to nearest lambda value.
        #
        # sci_scale $w.lc.lambda_corner_scale -from [set $this-lambda_min] -to [set $this-lambda_max] \
        #   -length [set $this-lambda_num] -resolution [set $this-lambda_resolution] -orient horizontal \
        #   -sliderlength 10 -bigincrement 0.5 -showvalue false -variable "$this-lambda_corner_from_scale" \
        #   -state disabled -foreground $color -command "$this update_lambda_corner"
        # pack $w.lc.lambda_corner_scale -side top -fill x
	
        global $this-lambda_lc
        sci_frame $w.lc.lam -relief flat
        sci_label $w.lc.lam.l -text "Lambda: "
        sci_entry $w.lc.lam.e -textvariable $this-lambda_lc -state disabled
        # sci_button $w.lc.lam.b -text "Update Lambda Value" -command "$this update_lambda_entry"
        pack $w.lc.lam.l $w.lc.lam.e -side left -padx 2 -pady 2 -fill x -expand yes
        pack $w.lc.lam -side top -fill x

        ######################
        #Execute Close Part
        ######################


        makeSciButtonPanel $w $w $this
 
        if {[set $this-reg_method] == "lcurve"} { $this uselcurve }
        if {[set $this-reg_method] == "single"} { $this entervalue }
        if {[set $this-reg_method] == "slider"} { $this useslider }

	moveToCursor $w
    }
}
