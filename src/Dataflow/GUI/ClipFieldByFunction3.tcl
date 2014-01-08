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

itcl::class SCIRun_NewField_ClipFieldByFunction3 {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name ClipFieldByFunction3
        set_defaults
    }

   method set_defaults {} {
      global $this-function
      global $this-help

      set $this-function "DATA1 < 0;"
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

        sci_labeledframe $w.inf -labeltext "Create Clipping Expression"
        set infoframe [$w.inf childsite]
        sci_frame $infoframe.info
        pack $infoframe.info -side left
        set info $infoframe.info
        sci_label $info.info1 -text "Function: expression(DATA,A,B,C,...)"
        sci_label $info.info2 -text "Input array: DATA1 (scalar/vector/tensor: data from field port) "
        sci_label $info.info2a -text "Input array: DATA2 (scalar/vector/tensor: data from field port) "
        sci_label $info.info2b -text "Input array: DATA3 (scalar/vector/tensor: data from field port) "
        sci_label $info.info3 -text "Input array: X1, Y1, Z1 (scalar: Cartensian coordinates of node/element)"
        sci_label $info.info3a -text "Input array: X2, Y2, Z2 (scalar: Cartensian coordinates of node/element)"
        sci_label $info.info3b -text "Input array: X3, Y3, Z3 (scalar: Cartensian coordinates of node/element)"
        sci_label $info.info4 -text "Input array: POS (vector: vector with node/element position)"
        sci_label $info.info5 -text "Input array: A, B, C, ... (scalar/vector/tensor: data from matrix ports)"
        sci_label $info.info6 -text "Input array: INDEX (scalar: number of the element)"
        sci_label $info.info7 -text "Input array: SIZE (scalar: number of elements)"
        sci_label $info.info8 -text "Input array: ELEMENT (element: object containing element properties)"

        grid $info.info1 -row 0 -column 0 -columnspan 2 -sticky w
        grid $info.info2 -row 1 -column 0 -sticky w
        grid $info.info2a -row 2 -column 0 -sticky w
        grid $info.info2b -row 3 -column 0 -sticky w
        grid $info.info3 -row 4 -column 0 -sticky w
        grid $info.info3a -row 5 -column 0 -sticky w
        grid $info.info3b -row 6 -column 0 -sticky w
        grid $info.info4 -row 1 -column 0 -sticky w
        grid $info.info5 -row 2 -column 1 -sticky w
        grid $info.info6 -row 3 -column 1 -sticky w
        grid $info.info7 -row 4 -column 1 -sticky w
        grid $info.info8 -row 5 -column 1 -sticky w

        pack $w.inf -side top -anchor w -fill x

        sci_labeledframe $w.ff -labeltext "Expression"
        set function [$w.ff childsite]
        option add *textBackground white	
        sci_scrolledtext $function.function -hscrollmode dynamic
        $function.function insert end [set $this-function]
        bind $function.function <Leave> "$this update_text"

        pack $w.ff -side top -anchor w -fill both -expand yes
        pack $function.function -side top -fill both -expand yes
  
        sci_labeledframe $w.clip -labeltext "Clipping location"
        set clip [$w.clip childsite]
        
        sci_radiobutton $clip.cell -text "Element Center" \
            -variable $this-method -value element
        sci_radiobutton $clip.nodeone -text "One Node" \
            -variable $this-method -value onenode
        sci_radiobutton $clip.nodemost -text "Most Nodes" \
            -variable $this-method -value mostnodes
        sci_radiobutton $clip.nodeall -text "All Nodes" \
            -variable $this-method -value allnodes
 
        pack $w.clip -side top -anchor w -fill x

        grid $clip.cell -row 0 -column 0 -sticky w
        grid $clip.nodeone -row 0 -column 1 -sticky w
        grid $clip.nodemost -row 0 -column 2 -sticky w
        grid $clip.nodeall -row 0 -column 3 -sticky w


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
    
}


