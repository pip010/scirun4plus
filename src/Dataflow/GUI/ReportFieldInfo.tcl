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


itcl::class SCIRun_MiscField_ReportFieldInfo {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name ReportFieldInfo

	# The width of the first column of the data display.
	setGlobal $this-firstwidth 12
    }

    method set_defaults {} {
	# These won't be saved.
      setGlobal $this-fldname "---"
      setGlobal $this-generation "---"
      setGlobal $this-typename "---"
      setGlobal $this-datamin "---"
      setGlobal $this-datamax "---"
      setGlobal $this-numnodes "---"
      setGlobal $this-numelems "---"
      setGlobal $this-dataat "---"
      setGlobal $this-cx "---"
      setGlobal $this-cy "---"
      setGlobal $this-cz "---"
      setGlobal $this-sizex "---"
      setGlobal $this-sizey "---"
      setGlobal $this-sizez "---"
      setGlobal $this-nodesx "---"
      setGlobal $this-nodesy "---"
      setGlobal $this-nodesz "---"
      setGlobal $this-geomsize "---"
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }
        sci_toplevel $w

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
        labelpair $att.l11 "Geometric size" $this-geomsize
        
        pack $att.l1 $att.l1a $att.l2 $att.l3 $att.l4 $att.l5 \
             $att.l7 $att.l8 $att.l9 $att.l10 $att.l11 -side top -expand yes -fill x

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




