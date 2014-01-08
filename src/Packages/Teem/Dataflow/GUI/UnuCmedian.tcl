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

#    File   : UnuCmedian.tcl
#    Author : Martin Cole
#    Date   : Mon Aug 25 10:14:23 2003

catch {rename Teem_UnuAtoM_UnuCmedian ""}

itcl::class Teem_UnuAtoM_UnuCmedian {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name UnuCmedian
        set_defaults
    }
    method set_defaults {} {
        global $this-radius
        global $this-weight
        global $this-bins
        global $this-pad
        

        set $this-radius 1
        set $this-weight 1.0
        set $this-bins 2048
        set $this-mode 0
        set $this-pad 0
    }

    method valid_int {new} {
        if {! [regexp "\\A\\d*\\Z" $new]} {
            return 0
        }
        return 1
          }

          method valid_float {new} {
        if {! [regexp "\\A\\d*\\.?\\d*\\Z" $new]} {
            return 0
        }
        return 1
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

        #radius
        sci_entryfield $w.f.options.radius -labeltext "Radius:" \
            -validate "$this valid_int %P"\
            -textvariable $this-radius
        #weight
        sci_entryfield $w.f.options.weight -labeltext "Weight:" \
            -validate "$this valid_float %P"\
            -textvariable $this-weight
        #bins
        sci_entryfield $w.f.options.bins -labeltext "Bins:" \
            -validate "$this valid_int %P"\
            -textvariable $this-bins
        #pad
        sci_label $w.f.options.padlabel -text "Pad:"
        sci_checkbutton $w.f.options.pad -variable $this-pad

        #pad
        sci_label $w.f.options.modelabel -text "Use Mode Filtering:"
        sci_checkbutton $w.f.options.mode -variable $this-mode	

        pack $w.f.options.radius $w.f.options.weight $w.f.options.bins -side top -expand yes -fill x
        pack $w.f.options.modelabel $w.f.options.mode -side left -anchor w 
        pack $w.f.options.padlabel $w.f.options.pad -side left -anchor w 
        

        pack $w.f -expand 1 -fill x
        makeSciButtonPanel $w $w $this
        moveToCursor $w

    }
}
