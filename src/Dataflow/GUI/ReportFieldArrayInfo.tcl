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


itcl::class SCIRun_FieldArray_ReportFieldArrayInfo {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name ReportArrayFieldInfo

      # The width of the first column of the data display.
      setGlobal $this-firstwidth 12
    }

    method update_range { } {
    
      set w .ui[modname]
      if {[winfo exists $w]} {
      
        upvar \#0 $this-selectfield_min min $this-selectfield_max max         
        $w.selectfield configure -from $min -to $max -variable $this-selectfield -command "$this-c needexecute"
      }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }
        sci_toplevel $w

        sci_scale $w.selectfield -from 0 -to 1 -variable $this-selectfield -showvalue true -orient horizontal -relief groove -label "Field Index:" 
        pack $w.selectfield -side top -anchor n -fill x -padx 8 -pady 5
        update_range

        sci_labeledframe $w.att -labelpos nw \
                         -labeltext "Input Field Attributes" 
                     
        pack $w.att -fill x -expand yes
        set att [$w.att childsite]
        
        entrypair $att.l1 "Name" $this-fldname
        entrypair $att.l1a "Generation" $this-generation
        labelpair $att.l2 "Type" $this-typename
        labelpair3 $att.l3 "Center (x,y,z)" $this-cx $this-cy $this-cz
        labelpair3 $att.l4 "Size (x,y,z)" $this-sizex $this-sizey $this-sizez
        labelpair2 $att.l5 "Data min,max" $this-datamin $this-datamax
        labelpair $att.l7 "# Nodes" $this-numnodes
        labelpair $att.l8 "# Elements" $this-numelems
        labelpair $att.l9 "Data at" $this-dataat
        labelpair3 $att.l10 "Dims (x,y,z)" $this-nodesx $this-nodesy $this-nodesz
        
        pack $att.l1 $att.l1a $att.l2 $att.l3 $att.l4 $att.l5 \
             $att.l7 $att.l8 $att.l9 $att.l10 -side top -expand yes -fill x

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }

    method entrypair { win text1 text2 } {

        sci_frame $win 
        pack $win -side top -padx 5
        sci_label $win.l1 -text $text1 -width [set $this-firstwidth] \
                -anchor w -just left
        sci_label $win.colon -text ":" -width 2 -anchor w -just left 

        sci_entry $win.l2 -textvar $text2 \
            -just left -width 40 \
            -relief flat -state disabled \
            -fore darkred -borderwidth 0 \
            -xscrollcommand [list $win.xscroll set]

        sci_scrollbar $win.xscroll -orient horizontal \
            -command [list $win.l2 xview]

        pack $win.l1 $win.colon $win.l2 -side left
        pack $win.xscroll -side left -fill x
    } 
}




