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


itcl::class SCIRun_Math_SelectSubMatrix {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name SelectSubMatrix
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        sci_checkbutton $w.rs -text "Select the following row interval" \
          -variable $this-row-select
        
        sci_frame $w.row
        sci_entryfield $w.row.start \
          -labeltext "Start index" \
          -textvariable $this-row-start -width 4       
        sci_entryfield $w.row.end \
          -labeltext "End index" \
          -textvariable $this-row-end -width 4       
        pack $w.row.start $w.row.end -side left
        
        pack $w.rs $w.row -anchor w

        sci_checkbutton $w.cs -text "Select the following column interval" \
          -variable $this-col-select
          
        sci_frame $w.col
        sci_entryfield $w.col.start \
          -labeltext "Start index" \
          -textvariable $this-col-start -width 4        
        sci_entryfield $w.col.end \
          -labeltext "End index" \
          -textvariable $this-col-end -width 4        
        pack $w.col.start $w.col.end -side left
        
        pack $w.cs $w.col -anchor w

        makeSciButtonPanel $w $w $this
        moveToCursor $w
     }
}
