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


itcl::class SCIRun_Visualization_GenerateStreamLines {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GenerateStreamLines
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w

        sci_checkbutton $w.autoparam -text "Automatically generate Error Tolerance and Step Size" -variable $this-auto-parameterize

        sci_frame $w.e -relief groove -borderwidth 2
        sci_frame $w.e.labs
        sci_frame $w.e.ents
        sci_label $w.e.labs.tolerance -text "Error Tolerance" -just left
        sci_entry $w.e.ents.tolerance -textvariable $this-tolerance
        sci_label $w.e.labs.stepsize -text "Step Size" -just left
        sci_entry $w.e.ents.stepsize -textvariable $this-stepsize
        sci_label $w.e.labs.maxsteps -text "Maximum Steps" -just left
        sci_entry $w.e.ents.maxsteps -textvariable $this-maxsteps
        sci_label $w.e.labs.nthreads -text "Number of Threads" -just left
        sci_spinint $w.e.ents.nthreads \
            -range {1 256} -step 1 \
            -textvariable $this-nthreads \
            -width 10 -fixed 10 -justify left


        pack $w.e.labs.tolerance $w.e.labs.stepsize $w.e.labs.maxsteps \
            $w.e.labs.nthreads -side top -anchor w
        pack $w.e.ents.tolerance $w.e.ents.stepsize $w.e.ents.maxsteps \
            $w.e.ents.nthreads -side top -anchor w

        pack $w.e.labs $w.e.ents -side left


        sci_frame $w.direction -relief groove -borderwidth 2
        sci_label $w.direction.label -text "Direction"
        sci_radiobutton $w.direction.neg -text "Negative" \
            -variable $this-direction -value 0
        sci_radiobutton $w.direction.both -text "Both" \
            -variable $this-direction -value 1
        sci_radiobutton $w.direction.pos -text "Positive" \
            -variable $this-direction -value 2

        pack $w.direction.label -side top -fill both
        pack $w.direction.neg $w.direction.both $w.direction.pos \
            -side left -fill both

        sci_frame $w.value -relief groove -borderwidth 2
        sci_label $w.value.label -text "Streamline Value"
        sci_frame $w.value.left
        sci_frame $w.value.middle
        sci_frame $w.value.right
        sci_radiobutton $w.value.left.value -text   "Seed Value" \
            -variable $this-value -value 0
        sci_radiobutton $w.value.left.index -text   "Seed Index" \
            -variable $this-value -value 1
        sci_radiobutton $w.value.middle.index -text "Integration Index" \
            -variable $this-value -value 2
        sci_radiobutton $w.value.middle.incr -text  "Integration Step" \
            -variable $this-value -value 3
        sci_radiobutton $w.value.right.delta -text  "Distance from Seed" \
            -variable $this-value -value 4
        sci_radiobutton $w.value.right.total -text  "Streamline Length" \
            -variable $this-value -value 5

	pack $w.value.left.value $w.value.left.index -side top -anchor w
	pack $w.value.middle.index $w.value.middle.incr -side top -anchor w
	pack $w.value.right.delta $w.value.right.total -side top -anchor w

        pack $w.value.left $w.value.middle $w.value.right \
	    -side left -fill both

        sci_frame $w.meth -relief groove -borderwidth 2
        sci_label $w.meth.label -text "Streamline Computation Method"
        sci_radiobutton $w.meth.cw -text "Cell Walk" \
            -variable $this-method -value 5
        sci_radiobutton $w.meth.ab -text "Adams-Bashforth Multi-Step" \
            -variable $this-method -value 0
        sci_radiobutton $w.meth.o2 -text "Adams Moulton Multi Step" \
            -variable $this-method \
            -value 1
        sci_radiobutton $w.meth.heun -text "Heun Method" \
            -variable $this-method -value 2
        sci_radiobutton $w.meth.rk4 -text "Classic 4th Order Runge-Kutta" \
            -variable $this-method -value 3
        sci_radiobutton $w.meth.rkf -text "Adaptive Runge-Kutta-Fehlberg" \
            -variable $this-method -value 4
             
        pack $w.meth.label -side top -fill both
        pack $w.meth.cw $w.meth.ab $w.meth.heun $w.meth.rk4 $w.meth.rkf \
         -side top -anchor w

        sci_checkbutton $w.filter -text "Filter Colinear Points" \
          -variable $this-remove-colinear-pts -justify left

        pack $w.meth $w.autoparam $w.e $w.direction $w.value $w.filter \
            -side top -e y -f both -padx 5 -pady 5

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}


