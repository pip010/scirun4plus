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
#    File   : FairMesh.tcl
#    Author : Martin Cole
#    Date   : Tue Mar 20 08:47:42 2007

itcl::class SCIRun_NewField_FairMesh {
    inherit Module
    
     constructor { {args ""} } {
        eval configure $args
        set name FairMesh
    }

    method ui {} {
	
        global $this-iterations
        global $this-method

        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        sci_frame $w.f     
        pack $w.f -expand yes -fill both

        make_labeled_radio $w.f.rb \
          "Weighting Method" "" top 1 \
          $this-method \
          { {"Fast (equal weights)" fast} {"Desbrun (curvature normal)" desbrun} }

        sci_entryfield $w.f.e -validate integer -width 8 \
          -labeltext "Iterations:" -textvariable $this-iterations

        sci_entryfield $w.f.l -validate real -width 8 \
          -labeltext "Relaxation Parameter:" -textvariable $this-lambda

        sci_entryfield $w.f.m -validate real -width 8 \
          -labeltext "Spatial cut off frequency:" -textvariable $this-mu


        pack $w.f.rb $w.f.e $w.f.l $w.f.m -anchor w -pady 2

        makeSciButtonPanel $w $w $this
        moveToCursor $w
     }

}


