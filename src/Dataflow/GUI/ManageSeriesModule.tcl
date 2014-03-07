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
 #  ManageSeriesModule.tcl: The ManageSeriesModule UI
 #  Written by:
 #   Allen R. Sanderson
 #   SCI Institute
 #   University of Utah
 #   March 2007
 ##

itcl::class ManageSeriesModule {
    inherit Module

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        ###### Save the num-ports because the iwidget resets it
        global $this-num-ports
        set quantity [set $this-num-ports]

        sci_spinint $w.ports -labeltext "Number of Ports to Manage: " \
            -range {1 4} -step 1 \
            -textvariable $this-num-ports \
            -width 10 -fixed 10 -justify right
        
        $w.ports delete 0 end
        $w.ports insert 0 $quantity

        sci_button $w.clear -text "Clear Ports" \
            -command "$this clear"

        pack $w.ports $w.clear -pady 5 -side top

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }

    method clear {} {
      $this-c clear
      $this-c needexecute
    }
}
