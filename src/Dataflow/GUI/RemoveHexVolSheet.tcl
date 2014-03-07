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
 #  RemoveHexVolSheet.tcl: The RemoveHexVolSheet UI
 #  Written by:
 #   Jason Shepherd
 #   Department of Computer Science
 #   University of Utah
 #   May 2006
 ##

catch {rename SCIRun_NewField_RemoveHexVolSheet ""}

itcl::class SCIRun_NewField_RemoveHexVolSheet {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name RemoveHexVolSheet
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return;
        }

        sci_toplevel $w
        wm minsize $w 150 80

        sci_frame $w.f1
        sci_label $w.f1.l -text "List of Edges:"
        sci_entry $w.f1.e -width 30 -text $this-edge-list
        bind $w.f1.e <Return> "$this-c needexecute"
        pack $w.f1.l $w.f1.e -side left -fill both -expand 1
        pack $w.f1 -fill x

        sci_frame $w.f
        sci_frame $w.fb
        pack $w.f $w.fb -padx 2 -pady 2 -side top -expand yes

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }
}
