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



catch {rename CoregisterPointClouds ""}

package require Iwidgets 3.0   

itcl::class SCIRun_MiscField_CoregisterPointClouds {
    inherit Module
    
     constructor { {args ""} } {
        eval configure $args
	set name CoregisterPointClouds

	global $this-method
	set $this-method "Analytic"
	trace variable $this-method w "$this switch_to_method"
    }
    
    method switch_to_method {name1 name2 op} {
      #puts stdout "switching"
      set window .ui[modname]
      if {[winfo exists $window]} {
          set mf [$window.f.meth childsite]
          $mf.tabs view [set $this-method]
      }
        }

        method ui {} {
      set w .ui[modname]
      if {[winfo exists $w]} {
          return
      }
      
      sci_toplevel $w
      sci_frame $w.f 
      pack $w.f -padx 2 -pady 2 -expand 1 -fill x
      set n "$this-c needexecute "
      
      sci_checkbutton $w.f.s -text "Allow scale" \
        -variable $this-allowScale
      sci_checkbutton $w.f.r -text "Allow rotate" \
        -variable $this-allowRotate
      sci_checkbutton $w.f.t -text "Allow translate" \
        -variable $this-allowTranslate
      
      pack $w.f.s $w.f.r $w.f.t

      #  Methods
      sci_labeledframe $w.f.meth -labelpos nw -labeltext "Methods"
      set mf [$w.f.meth childsite]
      
      sci_tabnotebook  $mf.tabs -raiseselect true 
      #-fill both
      pack $mf.tabs -side top

      #  Method:

      set alg [$mf.tabs add -label "Analytic" \
        -command "$this select-alg 0"]
      
      set alg [$mf.tabs add -label "Procrustes" \
        -command "$this select-alg 1"]

      set alg [$mf.tabs add -label "Simplex" \
        -command "$this select-alg 2"]

      sci_frame $alg.s
      pack $alg.s -side top
      sci_label $alg.s.l -text "Seed: "
      sci_entry $alg.s.e -textvariable $this-seed
      bind $alg.s.e <KeyPress-Return> $n
      pack $alg.s.l $alg.s.e -side left -fill x -expand 1

      sci_frame $alg.iters
      pack $alg.iters -side top
      sci_label $alg.iters.l -text "Max iterations: "
      sci_entry $alg.iters.e -textvariable $this-iters
      bind $alg.iters.e <KeyPress-Return> $n
      pack $alg.iters.l $alg.iters.e -side left -fill x -expand 1

      sci_frame $alg.misfit
      pack $alg.misfit -side top
      sci_label $alg.misfit.l -text "Misfit tolerance: "
      sci_entry $alg.misfit.e -textvariable $this-misfitTol
      bind $alg.misfit.e <KeyPress-Return> $n
      pack $alg.misfit.l $alg.misfit.e -side left -fill x -expand 1

      sci_button $alg.b -text "Abort" -command "$this-c stop"
      pack $alg.b -side top

      $mf.tabs view [set $this-method]
      $mf.tabs configure -tabpos "n"
      
      pack $mf.tabs -side top
      pack $w.f.meth -side top

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }

    method select-alg { alg } {
      global $this-method

      if { $alg == 0 } {
          set $this-method "Analytic"
      } elseif { $alg == 1 } {
          set $this-method "Procrustes"
      } else {
          set $this-method "Simplex"
      }
    }

}
