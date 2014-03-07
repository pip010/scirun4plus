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


    catch {rename BioPSE_Inverse_SolveInverseProblemWithTikhonovSVD ""}

    itcl::class BioPSE_Inverse_SolveInverseProblemWithTikhonovSVD {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name SolveInverseProblemWithTikhonovSVD
        set_defaults
    }
    method set_defaults {} {
        global $this-tex_var
        global $this-lambda_fix
        global $this-lambda_sld
        global $this-lambda_lc
        global $this-reg_method
        global $this-lambda_num
        global $this-have_ui
        global $this-lambda_min
        global $this-lambda_max
        set $this-tex_var 0.02
        set $this-lambda_fix 0.02
        set $this-lambda_sld 0.0
        set $this-lambda_lc 0.0
            set $this-reg_method lcurve
        set $this-lambda_num 200
        set $this-have_ui 0
        set $this-lambda_min 1e-6
        set $this-lambda_max 10
    }
    method entervalue {} {
        set w .ui[modname]
        set color "#505050"
        #$this-c needexecute
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
        $this calc_log [set $this-lambda_max] [set $this-lambda_min] [set $this-lambda_sld]
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
        #$this-c needexecute
        $w.sld.sld.f3.s configure -state disabled -foreground $color
        $w.ent.f2.l1 configure -foreground $color
        $w.ent.f2.e1 configure -state disabled -foreground $color
        $w.lc.lam.l configure -foreground black
        $w.lc.lam.e configure -state normal -foreground black
        $w.sld.sld.f2.lvalue configure -state disabled -foreground $color
        $w.sld.sld.f2.l1 configure -foreground $color
        $w.lc.f2.e1 configure -state normal -foreground black
        $w.lc.f2.l1 configure -foreground black
    }
    method execrunmode {} {
        set w .ui[modname]
        $this-c needexecute
    }
    ##############
    method plot_graph {lcurv lcor lam} {
        set w .ui[modname]
        if {![winfo exists $w]} {
            return
        }
        global $this-lambda_lc
        set $this-lambda_lc $lam

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
        $w.lc.data.g element delete LCurve
        $w.lc.data.g element delete LCorner
		$w.lc.data.g element create LCurve -data "1 0"
        $w.lc.data.g element create LCorner -data "1 0"
	}
    ##############
  	method calc_log {l_max l_min slider_value} {
        global $this-tex_var
        set w .ui[modname]
        set $this-tex_var [expr {[set $this-lambda_min]*pow(10,$slider_value/5*log10([set $this-lambda_max]/[set $this-lambda_min]))}] 
	}	

    ##############	
    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return;
        }

        global $this-have_ui
        set $this-have_ui 1

        ###################################
        #THIS PART IS FOR MAKING THE TITLE
        ###################################
        sci_toplevel $w
        wm minsize $w 150 20
        sci_frame $w.titlefr -relief groove -border 3.5
        sci_label $w.titlefr.l -text "Select Method for Lambda"
        pack $w.titlefr -side top -padx 2 -pady 2 -fill both
        pack $w.titlefr.l -side top -fill x

        global $this-lambda_fix
        global $this-lambda_sld
        global $this-reg_method
        
        #######################
        # Entry radio-button
        #######################
        sci_frame $w.ent -relief groove -border 3
        pack $w.ent -side top -expand yes -fill x

        sci_frame $w.ent.f1 -relief flat
        pack $w.ent.f1 -side top -expand yes -fill x

        sci_radiobutton $w.ent.f1.b -text "Enter value" \
          -variable "$this-reg_method" -value single \
          -command "$this entervalue"

        pack $w.ent.f1.b -side left

        sci_frame $w.ent.f2 -relief flat
        pack $w.ent.f2 -side top -expand yes -fill x
        sci_label $w.ent.f2.l1 -text "Lambda: "
        sci_entry $w.ent.f2.e1 -textvariable $this-lambda_fix
        pack $w.ent.f2.l1 $w.ent.f2.e1 -side left -expand yes \
          -fill x -padx 2 -pady 2
        #bind $w.ent.f2.e1 <Return> $n
        #########################
        # Slider radio-button
        #########################
        sci_frame $w.sld -relief flat
        pack $w.sld -side top -expand yes -fill x
	
        sci_frame $w.sld.sld -relief groove -border 3
        pack $w.sld.sld -side left -expand yes -fill x
        sci_frame $w.sld.rng -relief groove -border 3
        pack $w.sld.rng -side right  -fill y
        
        sci_frame $w.sld.rng.f1 -relief flat
        pack  $w.sld.rng.f1 -side top -fill x
        sci_label $w.sld.rng.f1.l1 -text "Lambda Range"	
        pack  $w.sld.rng.f1.l1 -side top -fill x

        sci_frame $w.sld.rng.f2 -relief flat
        pack  $w.sld.rng.f2 -side top -fill x
        sci_entry $w.sld.rng.f2.e1 -textvariable $this-lambda_min
        sci_label $w.sld.rng.f2.l1 -text "From "
        pack $w.sld.rng.f2.e1  $w.sld.rng.f2.l1 -side right
        
        sci_frame $w.sld.rng.f3 -relief flat
        pack  $w.sld.rng.f3 -side top -fill x
        sci_entry $w.sld.rng.f3.e1 -textvariable $this-lambda_max
        sci_label $w.sld.rng.f3.l1 -text "  To  "
        pack $w.sld.rng.f3.e1  $w.sld.rng.f3.l1 -side right

        sci_frame $w.sld.sld.f1 -relief flat
        pack $w.sld.sld.f1 -side top -expand yes -fill x
        sci_radiobutton $w.sld.sld.f1.b -text "Choose using slider" \
          -variable "$this-reg_method" -value slider \
          -command "$this useslider"
        pack $w.sld.sld.f1.b -side left

        sci_frame $w.sld.sld.f2 -relief flat
        pack $w.sld.sld.f2 -side top -expand yes -fill x
	
        sci_entry $w.sld.sld.f2.lvalue -textvariable "$this-tex_var" -state disabled       
        sci_label $w.sld.sld.f2.l1 -text "Lambda: "
        pack $w.sld.sld.f2.l1 $w.sld.sld.f2.lvalue -side left
         
        sci_frame $w.sld.sld.f3 -relief flat
              pack $w.sld.sld.f3 -side top -expand yes -fill x
        sci_scale $w.sld.sld.f3.s -from 0.0 -to 5.0 \
          -resolution 0.01 -orient horizontal \
          -variable "$this-lambda_sld" -showvalue false \
          -command "$this calc_log [set $this-lambda_max] [set $this-lambda_min] " 
        pack $w.sld.sld.f3.s -side left -expand yes -fill x -padx 2 -pady 2
	
        #######################
        # LCurve radio-button
        #######################

        sci_frame $w.lc -relief groove -border 3
        pack $w.lc -side top -expand yes -fill x
	
        sci_frame $w.lc.f1 -relief flat
        pack $w.lc.f1 -side top -expand yes -fill x
        sci_radiobutton $w.lc.f1.b -text "L-curve             " \
          -variable "$this-reg_method" -value lcurve \
          -command "$this uselcurve"
	
        pack $w.lc.f1.b -side left
	
        sci_frame $w.lc.f2 -relief flat
        pack $w.lc.f2 -side top -expand yes -fill x
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
        pack $w.lc.data.g -side top -fill x
        pack $w.lc.data -side top -fill x
		
				
        global $this-lambda_lc
        sci_frame $w.lc.lam -relief flat
        sci_label $w.lc.lam.l -text "Lambda Corner: "
        sci_entry $w.lc.lam.e -textvariable $this-lambda_lc -state disabled
        pack $w.lc.lam.l $w.lc.lam.e -side left -fill x -expand 1
        pack $w.lc.lam -side top -fill x

        ######################
        #Execute Close Part
        ######################
        makeSciButtonPanel $w $w $this
        moveToCursor $w
	
        if {[set $this-reg_method] == "lcurve"} { $this uselcurve }
        if {[set $this-reg_method] == "single"} { $this entervalue }
        if {[set $this-reg_method] == "slider"} { $this useslider }
    }
}
