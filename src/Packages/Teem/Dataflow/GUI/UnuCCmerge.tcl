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

#    File   : UnuCCmerge.tcl
#    Author : Darby Van Uitert
#    Date   : April 2004

itcl::class Teem_UnuAtoM_UnuCCmerge {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name UnuCCmerge
        set_defaults
    }

    method set_defaults {} {
        global $this-dir
        set $this-dir 0

        global $this-maxsize
        set $this-maxsize 0

        global $this-maxneigh
        set $this-maxneigh 1

        global $this-connectivity
        set $this-connectivity 1
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w

        sci_frame $w.f
        pack $w.f -padx 2 -pady 2 -side top -expand yes
        
        sci_frame $w.f.options
        pack $w.f.options -side top -expand yes

        sci_labeledframe $w.f.options.dir \
            -labeltext "Value Driven Merging" \
            -labelpos nw
        pack $w.f.options.dir -side top -expand yes -fill x

        set dir [$w.f.options.dir childsite]

        sci_radiobutton $dir.dir1 \
            -text "Default - Merging can go either way" \
            -variable $this-dir \
            -value {0}
        
        sci_radiobutton $dir.dir2 \
            -text "Dark islands get merged with bright surrounds" \
            -variable $this-dir \
            -value 1

        sci_radiobutton $dir.dir3 \
            -text "Bright surrounds get merged with dark islands" \
            -variable $this-dir \
            -value {-1}

        pack $dir.dir1 $dir.dir2 $dir.dir3 -side top -anchor nw
        sci_entryfield $w.f.options.maxsize \
          -labeltext "Max Size:" -textvariable $this-maxsize
        pack $w.f.options.maxsize -side top -expand yes -fill x


        sci_entryfield $w.f.options.maxneigh \
          -labeltext "Max Neighbors:" -textvariable $this-maxneigh
        pack $w.f.options.maxneigh -side top -expand yes -fill x

        sci_entryfield $w.f.options.connectivity \
          -labeltext "Connectivity:" -textvariable $this-connectivity
        pack $w.f.options.connectivity  -side top -expand yes -fill x


        pack $w.f -expand 1 -fill x
        makeSciButtonPanel $w $w $this
        moveToCursor $w

    }
}
