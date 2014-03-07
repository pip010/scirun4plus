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

#    File   : UnuSort.tcl
#    Author : Allen R. Sanderson
#    Date   : March 2007

catch {rename Teem_UnuNtoZ_UnuSort ""}

itcl::class Teem_UnuNtoZ_UnuSort {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name UnuSort
        set_defaults
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        global $this-index
        global $this-min
        global $this-max

        set index_tmp [set $this-index]
        set min_tmp [set $this-min]
        set max_tmp [set $this-max]

        sci_spinint $w.index -labeltext "Index for sorting: " \
            -range "$min_tmp $max_tmp" -step 1 \
            -textvariable $this-index \
            -width 10 -fixed 10 -justify right
        
        pack $w.index -side top -expand 1 -fill x -padx 5

        $w.index delete 0 end
        $w.index insert 0 $index_tmp

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }

    method set_min_max { } {
        set w .ui[modname]

        global $this-index
        global $this-min
        global $this-max

        set index_tmp [set $this-index]
        set min_tmp [set $this-min]
        set max_tmp [set $this-max]

        if [ expr [winfo exists $w] ] {

                  if { $max_tmp == $min_tmp } {
                      set max_tmp [expr $min_tmp + 1]
                  }

            $w.index configure -range "$min_tmp $max_tmp"
            $w.index delete 0 end
            $w.index insert 0 $index_tmp
        }
    }
}
