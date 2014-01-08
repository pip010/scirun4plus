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

itcl::class SCIRun_DataArrayMath_CalculateDataArray {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name CalculateDataArray
        set_defaults
    }

    method set_defaults {} {
      global $this-function
      global $this-help
      global $this-format

      set $this-function "RESULT = abs(DATA);"
      set $this-help ""
      set $this-format "Scalar"
      
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

        sci_labeledframe $w.inf -labeltext "Create Data Array"
        set infoframe [$w.inf childsite]
        sci_frame $infoframe.info
        pack $infoframe.info -side left
        set info $infoframe.info
        sci_label $info.info1 -text "Function: RESULT = function(DATA,A,B,C,...)"
        sci_label $info.info2 -text "Input array: DATA, A, B, C, ... (scalar,vector, or tensor)"
        sci_label $info.info3 -text "Output array: RESULT (scalar)"
        sci_label $info.info4 -text "Element index: INDEX (scalar)"
        sci_label $info.info5 -text "Number of elements: SIZE (scalar)"
        grid $info.info1 -row 0 -column 0 -sticky w
        grid $info.info2 -row 1 -column 0 -sticky w
        grid $info.info3 -row 2 -column 0 -sticky w
        grid $info.info4 -row 3 -column 0 -sticky w
        grid $info.info5 -row 4 -column 0 -sticky w
        pack $w.inf -side top -anchor w -fill x

        sci_labeledframe $w.of -labeltext "output type"
        set otype [$w.of childsite]
        pack $w.of -side top -anchor w -fill x
        
        labelcombo $otype.otype "Field Output Data Type" \
          {Scalar Vector Tensor} \
          $this-format
        
        sci_labeledframe $w.ff -labeltext "function"
        set function [$w.ff childsite]
        option add *textBackground white	
        sci_scrolledtext $function.function -height 60 -hscrollmode dynamic
        $function.function insert end [set $this-function]
        bind $function.function <Leave> "$this update_text"

        pack $w.ff -side top -anchor w -fill both 
        pack $function.function -side top -fill both 

        sci_button $w.help -text "Available Functions" -command "$this showhelp"
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

      sci_labeledframe $w.hf -labeltext "available functions"
      set help [$w.hf childsite]
      option add *textBackground white	
      iwidgets::scrolledhtml $help.help -height 60 -hscrollmode dynamic -width 500p -height 300p        
      $help.help render [set $this-help]
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
      for {set elem [lindex $arglist $i]} {$i<$length} {incr i 1; set elem [lindex $arglist $i]} {
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


