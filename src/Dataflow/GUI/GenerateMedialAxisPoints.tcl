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


itcl::class SCIRun_NewField_GenerateMedialAxisPoints {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GenerateMedialAxisPoints
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }
        sci_toplevel $w

        global $this-refinement-levels
        sci_entryfield $w.levels -labeltext "Number of refinement levels =" -textvariable $this-refinement-levels

        global $this-axis-min-angle
        sci_entryfield $w.axisangle -labeltext "Minimum angle between medial axis rays =" -textvariable $this-axis-min-angle

        global $this-normal-min-angle
        sci_entryfield $w.normalangle -labeltext "Minimum angle between surfaces =" -textvariable $this-normal-min-angle

        global $this-max-distance
        sci_entryfield $w.maxdistance -labeltext "Maximum distance difference to qualify for medial axis =" -textvariable $this-max-distance

        global $this-forward-multiplier
        sci_entryfield $w.forwardmult -labeltext "Value of forward projection to surface multiplier=" -textvariable $this-forward-multiplier

        global $this-backward-multiplier
        sci_entryfield $w.backmult -labeltext "Value of backward projection to surface multiplier=" -textvariable $this-backward-multiplier

        pack $w.levels $w.axisangle $w.normalangle $w.maxdistance $w.forwardmult $w.backmult -side top -expand yes -fill x

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
