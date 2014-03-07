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

catch {rename SCIRun_ChangeMesh_ResampleRegularMesh ""}

itcl::class SCIRun_ChangeMesh_ResampleRegularMesh {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name ResampleRegularMesh
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        sci_labeledframe $w.dims -labeltext "Dimensions of output mesh"
        pack $w.dims -fill x -expand yes -side top
        set dims [$w.dims childsite]
        
        sci_label $dims.xlabel -text "X-axis"
        sci_label $dims.ylabel -text "Y-axis"
        sci_label $dims.zlabel -text "Z-axis"
        
        sci_entry $dims.xentry -textvariable $this-xdim
        sci_entry $dims.yentry -textvariable $this-ydim
        sci_entry $dims.zentry -textvariable $this-zdim
        
        sci_label $dims.note -text "Enter number of samples or scaling factor (e.g. `x0.5')"
        grid $dims.xlabel -column 0 -row 0 -sticky w
        grid $dims.ylabel -column 0 -row 1 -sticky w
        grid $dims.zlabel -column 0 -row 2 -sticky w

        grid $dims.xentry -column 1 -row 0 -sticky w
        grid $dims.yentry -column 1 -row 1 -sticky w
        grid $dims.zentry -column 1 -row 2 -sticky w
        
        grid $dims.note -columnspan 2 -column 0 -row 3 -sticky w
        
        sci_labeledframe $w.kernel -labeltext "Resampling Kernel"
        pack $w.kernel -fill x -expand yes -side top
        set kernel [$w.kernel childsite]  
  
        sci_radiobutton $kernel.b1 -text "Box" -value "box" -variable $this-method
        sci_radiobutton $kernel.b2 -text "Tent" -value "tent" -variable $this-method
        sci_radiobutton $kernel.b3 -text "Cubic (Catmull-Rom)" -value "cubiccr" -variable $this-method
        sci_radiobutton $kernel.b4 -text "Cubic (B-Spline)" -value "cubicbs" -variable $this-method
        sci_radiobutton $kernel.b5 -text "Quartic" -value "quartic" -variable $this-method
        sci_radiobutton $kernel.b6 -text "Gaussian" -value "gaussian" -variable $this-method

        grid $kernel.b1 -row 0 -column 0 -sticky w
        grid $kernel.b2 -row 1 -column 0 -sticky w
        grid $kernel.b3 -row 2 -column 0 -sticky w
        grid $kernel.b4 -row 0 -column 1 -sticky w
        grid $kernel.b5 -row 1 -column 1 -sticky w
        grid $kernel.b6 -row 2 -column 1 -sticky w

        sci_label $kernel.lsigma -text "Gaussian sigma"
        sci_entry $kernel.sigma -textvariable $this-sigma

        sci_label $kernel.lextend -text "Gaussian extend"
        sci_entry $kernel.extend -textvariable $this-extend

        grid $kernel.lsigma -row 1 -column 2
        grid $kernel.sigma -row 1 -column 3

        grid $kernel.lextend -row 2 -column 2
        grid $kernel.extend -row 2 -column 3

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
