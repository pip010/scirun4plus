#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  License for the specific language governing rights and limitations under
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


# GUI for FiniteTimeLyaponovExponent module
# by Allen R. Sanderson
# SCI Institute
# University of Utah
# September 2005

catch {rename FusionPSE_Fields_FiniteTimeLyaponovExponent ""}

itcl::class FusionPSE_Fields_FiniteTimeLyaponovExponent {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name FiniteTimeLyaponovExponent
    }
    
    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }

	toplevel $w

	frame $w.stepsize
	label $w.stepsize.label -text "Step size:"
	entry $w.stepsize.entry -width 20 -text $this-stepsize

	pack $w.stepsize.label $w.stepsize.entry \
	    -side left -anchor w


	frame $w.numsteps
	label $w.numsteps.label -text "Num of Steps:"
	entry $w.numsteps.entry -width 20 -text $this-numsteps

	pack $w.numsteps.label $w.numsteps.entry \
	    -side left -anchor w


	frame $w.delta
	label $w.delta.label -text "Delta:"
	entry $w.delta.entry -width 20 -text $this-delta

	pack $w.delta.label $w.delta.entry \
	    -side left -anchor w


	frame $w.dimension -relief groove -borderwidth 2
	label $w.dimension.label -text "Dimension"

	radiobutton $w.dimension.two     -text "2D" \
	    -variable $this-dim -value 2
	radiobutton $w.dimension.three   -text "3D" \
	    -variable $this-dim -value 3

	pack $w.dimension.label -side top -fill both
	pack $w.dimension.two $w.dimension.three -side left -anchor w


	frame $w.plane -relief groove -borderwidth 2
	label $w.plane.label -text "Plane"

	radiobutton $w.plane.xy     -text "X-Y" \
	    -variable $this-plane -value 2
	radiobutton $w.plane.xz   -text "X-Z" \
	    -variable $this-plane -value 1
	radiobutton $w.plane.yz   -text "Y-Z" \
	    -variable $this-plane -value 0

	pack $w.plane.label -side top -fill both
	pack $w.plane.xy $w.plane.xz $w.plane.yz -side left -anchor w


	pack $w.stepsize $w.numsteps $w.delta $w.dimension $w.plane \
	    -side top -fill x

	makeSciButtonPanel $w $w $this
	moveToCursor $w
    }
}
