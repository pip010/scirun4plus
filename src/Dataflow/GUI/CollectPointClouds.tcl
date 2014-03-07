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
 #  CollectPointClouds.tcl: The CollectPointClouds UI
 #
 #  Written by:
 #   Allen R. Sanderson
 #   SCI Institute
 #   University of Utah
 #   July 2007
 ##

catch {rename SCIRun_MiscField_CollectPointClouds ""}

itcl::class SCIRun_MiscField_CollectPointClouds {

    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name CollectPointClouds
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        ###### Save the num_fields because the iwidget resets it
        global $this-num_fields
        set quantity [set $this-num_fields]

        sci_spinint $w.values \
            -labeltext "Number of Fields to Collect: " \
            -range {2 999} -step 1 \
            -textvariable $this-num_fields \
            -width 10 -fixed 10 -justify right
        
        $w.values delete 0 end
        $w.values insert 0 $quantity


        sci_frame $w.misc

        sci_label $w.misc.label -text "Count:"
        sci_entry $w.misc.count -width 5 -textvariable $this-count -state disabled
        
        sci_button $w.misc.clear -text " Clear Data " \
            -command "$this-c clear"

        pack $w.misc.label $w.misc.count $w.misc.clear -padx 5 -side left

        pack $w.values $w.misc -pady 5 -side top

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
