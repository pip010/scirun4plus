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

itcl::class SCIRun_ChangeFieldData_CalculateFieldDataMetric {
   inherit Module

     constructor { {args ""} } {
        eval configure $args
      set name CalculateFieldDataMetric
  
      # The width of the first column of the data display.
      setGlobal $this-firstwidth 12
    }

    method update_text {} {
      set w .ui[modname]
      if {[winfo exists $w]} {
        set transform [$w.transform childsite]
        set $this-function [$transform.t.function get 1.0 end]
        }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        sci_labeledframe $w.transform -labeltext "Transform FieldData"
        set transform [$w.transform childsite]
        pack $w.transform -side top -anchor w -fill x     

        sci_frame $transform.t
        pack $transform.t -side top -anchor w  -fill x    

        sci_checkbutton $transform.t.enable -variable $this-enable-function -text "Transform the data before computing metric"
        grid $transform.t.enable -row 0 -column 0 -columnspan 2 -sticky w -pady 5    

        
        sci_label $transform.t.info1 -text "Function: RESULT = function(DATA,POS,ELEMENT,INDEX,...)"
        sci_label $transform.t.info2 -text "Input array: DATA (scalar/vector/tensor: data from field port) "
        sci_label $transform.t.info3 -text "Input array: X, Y, Z (scalar: Cartensian coordinates of node/element)"
        sci_label $transform.t.info4 -text "Input array: POS (vector: vector with node/element position)"
        sci_label $transform.t.info5 -text "Output array: RESULT (scalar/vector/tensor)"

        grid $transform.t.info1 -row 1 -column 0 -columnspan 2 -sticky w
        grid $transform.t.info2 -row 2 -column 0 -columnspan 2 -sticky w
        grid $transform.t.info3 -row 3 -column 0 -columnspan 2 -sticky w
        grid $transform.t.info4 -row 4 -column 0 -columnspan 2 -sticky w
        grid $transform.t.info5 -row 5 -column 0 -columnspan 2 -sticky w

        sci_label $transform.t.expression -text "Expression"
        grid $transform.t.expression -row 6 -column 0 -sticky w -pady 5

        sci_button $transform.t.help -text "Parser Help" -command "$this showhelp"
        grid $transform.t.help -row 6 -column 1 -sticky e   

        sci_scrolledtext $transform.t.function -height 60 -hscrollmode dynamic
        $transform.t.function insert end [set $this-function]
        bind $transform.t.function <Leave> "$this update_text"

        grid $transform.t.function -row 7 -column 0 -columnspan 2 -sticky we

        sci_labeledframe $w.metric -labeltext "Metric"
        set metric [$w.metric childsite]
        pack $w.metric -side top -anchor e -fill x
          
        sci_radiobutton $metric.min -text "Minimum" \
            -variable $this-method -value "min"
        grid $metric.min -row 0 -column 0 -sticky w

        sci_radiobutton $metric.max -text "Maximum" \
            -variable $this-method -value "max"
        grid $metric.max -row 0 -column 1 -sticky w

        sci_radiobutton $metric.median -text "Median" \
            -variable $this-method -value "median"
        grid $metric.median -row 0 -column 2 -sticky w

        sci_radiobutton $metric.valuemean -text "Value-Mean" \
            -variable $this-method -value "value-mean"
        grid $metric.valuemean -row 0 -column 3 -sticky w

        sci_radiobutton $metric.geommean -text "Geom-Mean" \
            -variable $this-method -value "geom-mean"
        grid $metric.geommean -row 1 -column 0 -sticky w
        
        sci_radiobutton $metric.sum -text "Sum" \
            -variable $this-method -value "sum"
        grid $metric.sum -row 1 -column 1 -sticky w
        
        sci_radiobutton $metric.integral -text "Integral" \
            -variable $this-method -value "integral"
        grid $metric.integral -row 1 -column 2 -sticky w
        
        sci_labeledframe $w.result -labeltext "Result"
        set result [$w.result childsite]
        pack $w.result -side top -anchor e -fill x
        
        labelpair4 $result.metric "Metric" $this-metric
        grid $result.metric -row 0 -column 0 -sticky w -pady 5 
        
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
    
}


