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


# GUI for GetFieldZeroCrossings module
# by Allen Sanderson
# April 2005

# This GUI interface is for selecting an axis and index for sub sampling a
# topologically structured field

catch {rename FusionPSE_Fields_GetFieldZeroCrossings ""}

itcl::class FusionPSE_Fields_GetFieldZeroCrossings {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GetFieldZeroCrossings
    }

    method ui {} {

	global $this-direction

        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }

        toplevel $w

	frame $w.axis

	label $w.axis.label -text "Axis" -width 6 -anchor w -just left

	radiobutton $w.axis.x -text "X" -width 3 \
		-anchor w -just left -variable $this-axis -value 0

	radiobutton $w.axis.y -text "Y" -width 3 \
		-anchor w -just left -variable $this-axis -value 1

	radiobutton $w.axis.z -text "Z" -width 3 \
		-anchor w -just left -variable $this-axis -value 2

	pack $w.axis.label $w.axis.x $w.axis.y \
	    $w.axis.z -side left


	pack $w.axis -side top
	
	makeSciButtonPanel $w $w $this
	moveToCursor $w
    }
}
