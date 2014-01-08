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


##
 #  BuildFEMatrix.tcl
 #  Written by:
 #   David Weinstein
 #   Department of Computer Science
 #   University of Utah
 #   Aug 1996, March 2001
 #
 # GUI code migrated from BioPSE_Forward_BuildFEMatrix
 ##

itcl::class SCIRun_FiniteElements_BuildFEMatrix {
    inherit Module
     constructor { {args ""} } {
       eval configure $args
       set name BuildFEMatrix
    }

    method ui {} {
      set w .ui[modname]

      if {[winfo exists $w]} {
        return
      }

      sci_toplevel $w

      sci_frame $w.np
      pack $w.np -side top -anchor nw -padx 4 -pady 2

      sci_label $w.np.l -text "Number of Threads"
      sci_entry $w.np.e -width 5 -textvariable $this-num-processors -justify center
      pack $w.np.l $w.np.e -side left -anchor w

      sci_checkbutton $w.b1 -text "Use Conductivity Basis Matrices" \
        -variable $this-use-basis
      sci_checkbutton $w.b2 -text "Force Symmetric Matrix" \
        -variable $this-force-symmetry
                        
      pack $w.b1 $w.b2 -side top -anchor nw -padx 4 -pady 2

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }
}
