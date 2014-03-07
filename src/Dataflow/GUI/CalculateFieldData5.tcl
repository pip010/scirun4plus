##
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

itcl::class SCIRun_ChangeFieldData_CalculateFieldData5 {
   inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name CalculateFieldData5
        set_defaults
    }

    method set_defaults {} {
      global $this-function
      global $this-format
      global $this-help

      set $this-function "RESULT = abs(DATA1);"
      set $this-format "Scalar"
      set $this-help ""
      
    }

    method update_text {} {
      set w .ui[modname]
      if {[winfo exists $w]} {
        set function [$w.ff childsite]
        set $this-function [$function.function get 1.0 end]
        }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        sci_labeledframe $w.inf -labeltext "Create New Field Data"
        set infoframe [$w.inf childsite]
        sci_frame $infoframe.info
        pack $infoframe.info -side left
        set info $infoframe.info
        sci_label $info.info1 -text "Function: RESULT = function(DATA,A,B,C,...)"
        sci_label $info.info2  -text "Input array: DATA1 (scalar/vector/tensor: data from field port) "
        sci_label $info.info2a -text "Input array: DATA2 (scalar/vector/tensor: data from field port) "
        sci_label $info.info2b -text "Input array: DATA3 (scalar/vector/tensor: data from field port) "
        sci_label $info.info2c -text "Input array: DATA4 (scalar/vector/tensor: data from field port) "
        sci_label $info.info2d -text "Input array: DATA5 (scalar/vector/tensor: data from field port) "
        sci_label $info.info3 -text "Input array: X, Y, Z (scalar: Cartensian coordinates of node/element)"
        sci_label $info.info4 -text "Input array: POS (vector: vector with node/element position)"
        sci_label $info.info5 -text "Input array: A, B, C, ... (scalar/vector/tensor: data from field data ports)"
        sci_label $info.info6 -text "Input array: INDEX (scalar: number of the element)"
        sci_label $info.info7 -text "Input array: SIZE (scalar: number of elements)"
        sci_label $info.info8 -text "Input array: ELEMENT (element: object containing element properties)"
        sci_label $info.info9 -text "Output array: RESULT (scalar)"

        grid $info.info1 -row 0 -column 0 -columnspan 2 -sticky w
        grid $info.info2 -row 1 -column 0 -sticky w
        grid $info.info2a -row 2 -column 0 -sticky w        
        grid $info.info2b -row 3 -column 0 -sticky w
        grid $info.info2c -row 4 -column 0 -sticky w
        grid $info.info2d -row 5 -column 0 -sticky w
        grid $info.info3 -row 6 -column 0 -sticky w
        grid $info.info4 -row 1 -column 1 -sticky w
        grid $info.info5 -row 2 -column 1 -sticky w
        grid $info.info6 -row 3 -column 1 -sticky w
        grid $info.info7 -row 4 -column 1 -sticky w
        grid $info.info8 -row 5 -column 1 -sticky w
        grid $info.info9 -row 6 -column 1 -sticky w

        pack $w.inf -side top -anchor w -fill x

        sci_frame $w.f
        pack $w.f -side top -anchor w -fill x

        sci_labeledframe $w.f.of -labeltext "Output Type"
        set otype [$w.f.of childsite]
        
        labelcombo $otype.otype "Output Data Type" \
          {Scalar Vector Tensor "Same as Input" "char" "unsigned char" \
					"short" "unsigned short" "int" "unsigned int" "float" "double"} \
          $this-format

        sci_labeledframe $w.f.cf -labeltext "Caching"
        set caching [$w.f.cf childsite]

        sci_checkbutton $caching.cache -text "Cache Result" \
          -variable $this-cache
        sci_button $caching.clear \
            -text "Clear Cache" -command "$this-c clear"

        sci_label $caching.label -text "Count:"

        sci_entry $caching.count -width 5 -textvariable $this-count  -state disabled

        pack $caching.cache -side left -padx 5 -pady 5
        pack $caching.label -side left -pady 5 -pady 5
        pack $caching.count -side left -padx 5 -pady 5
        pack $caching.clear -side left -padx 5 -pady 5
        pack $w.f.cf $w.f.of -side left -fill x -anchor w

        sci_labeledframe $w.ff -labeltext "function"
        set function [$w.ff childsite]
        option add *textBackground white	
        sci_scrolledtext $function.function -hscrollmode dynamic -height 30
        $function.function insert end [set $this-function]
        bind $function.function <Leave> "$this update_text"
                
        pack $w.ff -side top -anchor w -fill both -expand yes
        pack $function.function -side top -fill both -expand yes

        sci_button $w.help -text "Parser Help" -command "$this showhelp"
        pack $w.help -side top -anchor e 

        makeSciButtonPanel $w $w $this
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


