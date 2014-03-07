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


itcl::class SCIRun_ChangeFieldData_QueryFieldData {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name QueryFieldData
	
        # Trace variable for optionmenu so that it will display
        # the correct value when opening a saved network.
        global $this-outputdatatype
        trace variable $this-outputdatatype w \
            "$this set_combobox .otype.c $this-outputdatatype"
          }

    method update_text {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
          set $this-function [$w.function.row1 get 1.0 end]
        }
    }

    method set_text {} {
	set w .ui[modname]
        if {[winfo exists $w]} {
	    $w.functio.row1 clear
	    $w.function.row1 insert end [set $this-function]
        }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }
        sci_toplevel $w

      # Output Data Type
        sci_frame $w.output -relief groove -borderwidth 2

        labelcombo $w.output.otype "Output Data Type" \
            {"port 0 input" \
           "unsigned char" "unsigned short" "unsigned int" \
           char short int float double Vector Tensor} \
            $this-outputdatatype

        pack $w.output.otype -side top -anchor w

        pack $w.output -side top -fill x -padx 5 -pady 5

      # Function
        sci_frame $w.function -borderwidth 2

        sci_label $w.function.info -text "F(x, y, z, v0, v1, ..., result, count):"


        sci_label $w.function.info1 -text "where 'x', 'y', and 'z' are the coordinate values, 'v0', 'v1', 'v2', etc."
        sci_label $w.function.info2 -text "are data values that correspond to each valid input field, and optionally"
        sci_label $w.function.info3 -text "'result' corresponds to the previous output field (if queried - otherwise zero) and"
        sci_label $w.function.info4 -text "'count' corresponds to the number of times the queried output field has been used."

        option add *textBackground white
        sci_scrolledtext $w.function.row1 -height 60 -hscrollmode dynamic

        $w.function.row1 insert end [set $this-function]

        pack $w.function.info  -side top -anchor w -padx 5 -pady 5
        pack $w.function.row1  -side top -fill both -expand true -padx 5
        pack $w.function.info1 -side top -anchor w -padx 5
        pack $w.function.info2 -side top -anchor w -padx 5
        pack $w.function.info3 -side top -anchor w -padx 5
        pack $w.function.info4 -side top -anchor w -padx 5

        pack $w.function -side top -fill both -expand true -padx 5 -pady 5

      # querying
        sci_labeledframe $w.querying -labelpos nw -labeltext "Querying"
        set querying [$w.querying childsite]

        sci_button $querying.clear \
            -text " Clear queried result " -command "$this-c clear"

        sci_label $querying.clabel -text "Count:"
        sci_label $querying.nlabel -text "Datasets:"

        sci_entry $querying.count -width 5 -textvariable $this-count  -state disabled
        sci_entry $querying.number_of_datasets -width 5 -textvariable $this-number_of_datasets  -state disabled
        
        pack $querying.clabel             -side left -pady  5
        pack $querying.count              -side left -padx 15 -pady 5
        pack $querying.nlabel             -side left -pady  5
        pack $querying.number_of_datasets -side left -padx 15 -pady 5
        pack $querying.clear              -side left -padx  5 -pady 5

        pack $w.querying -side top -fill x -padx 5 -pady 5
 
        sci_button $w.help -text "Parser Help" -command "$this showhelp"
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

    method labelcombo { win text1 arglist var} {
      sci_frame $win 
      pack $win -side top -padx 5
      sci_label $win.l1 -text $text1 \
              -anchor w -just left
      sci_label $win.colon  -text ":" -width 2 -anchor w -just left
      sci_optionmenu $win.c -foreground darkred \
        -command " $this comboget $win.c $var "

      set i 0
      set found 0
      set length [llength $arglist]
      for {set elem [lindex $arglist $i]} {$i<$length} \
          {incr i 1; set elem [lindex $arglist $i]} {
          if {"$elem"=="[set $var]"} {
        set found 1
          }
          $win.c insert end $elem
      }

      if {!$found} {
          $win.c insert end [set $var]
      }

      $win.c select [set $var]

      sci_label $win.l2 -text "" -width 20 -anchor w -just left

      # hack to associate optionmenus with a textvariable
      # bind $win.c <Map> "$win.c select {[set $var]}"

      pack $win.l1 $win.colon -side left
      pack $win.c $win.l2 -side left	
    }

    method comboget { win var } {
      if {![winfo exists $win]} {
          return
      }
      if { "$var"!="[$win get]" } {
          set $var [$win get]
      }
    }

    method set_combobox { win var name1 name2 op } {
      set w .ui[modname]
      set menu $w.$win
      if {[winfo exists $menu]} {
          $menu select $var
      }
    }
}


