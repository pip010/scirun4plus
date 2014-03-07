##
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


itcl::class SCIRun_Bundle_InsertColorMapsIntoBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name InsertColorMapsIntoBundle
    }

    method ui {} {
        
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        # input matrix names

        global $this-colormap1-name
        global $this-colormap2-name
        global $this-colormap3-name
        global $this-colormap4-name
        global $this-colormap5-name
        global $this-colormap6-name
        global $this-bundlename

        sci_toplevel $w 

        wm minsize $w 100 150

        sci_labeledframe $w.frame -labeltext "COLORMAP INPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill x
        sci_frame $w.frame2
        pack $w.frame2 -fill x
        sci_label $w.frame2.label -text "Name bundle object :"
        sci_entry $w.frame2.entry -textvariable $this-bundlename
        pack $w.frame2.label -side left 
        pack $w.frame2.entry -side left -fill x

        sci_tabnotebook $childframe.pw -height 100 -width 450 -tabpos n
        $childframe.pw add -label "ColorMap1"
        $childframe.pw add -label "ColorMap2" 
        $childframe.pw add -label "ColorMap3" 
        $childframe.pw add -label "ColorMap4"
        $childframe.pw add -label "ColorMap5" 
        $childframe.pw add -label "ColorMap6" 
        $childframe.pw select 0

        pack $childframe.pw -fill x -expand yes

        set colormap1 [$childframe.pw childsite 0]
        set colormap2 [$childframe.pw childsite 1]
        set colormap3 [$childframe.pw childsite 2]
        set colormap4 [$childframe.pw childsite 3]
        set colormap5 [$childframe.pw childsite 4]
        set colormap6 [$childframe.pw childsite 5]

        sci_frame $colormap1.name
        sci_frame $colormap1.options
        pack $colormap1.name $colormap1.options -side top -fill x -expand yes -padx 5p

        sci_label $colormap1.name.label -text "Name"
        sci_entry $colormap1.name.entry -textvariable $this-colormap1-name
        pack $colormap1.name.label -side left 
        pack $colormap1.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $colormap1.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $colormap1.options.replace -side left

        sci_frame $colormap2.name
        sci_frame $colormap2.options
        pack $colormap2.name $colormap2.options -side top -fill x -expand yes -padx 5p

        sci_label $colormap2.name.label -text "Name"
        sci_entry $colormap2.name.entry -textvariable $this-colormap2-name
        pack $colormap2.name.label -side left 
        pack $colormap2.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $colormap2.options.replace -variable $this-replace2 -text "Replace field, if it already exists"
        pack $colormap2.options.replace -side left

        sci_frame $colormap3.name
        sci_frame $colormap3.options
        pack $colormap3.name $colormap3.options -side top -fill x -expand yes -padx 5p

        sci_label $colormap3.name.label -text "Name"
        sci_entry $colormap3.name.entry -textvariable $this-colormap3-name
        pack $colormap3.name.label -side left 
        pack $colormap3.name.entry -side left -fill x -expand yes

        sci_checkbutton $colormap3.options.replace -variable $this-replace3 -text "Replace field, if it already exists"
        pack $colormap3.options.replace -side left

        sci_frame $colormap4.name
        sci_frame $colormap4.options
        pack $colormap4.name $colormap4.options -side top -fill x -expand yes -padx 5p

        sci_label $colormap4.name.label -text "Name"
        sci_entry $colormap4.name.entry -textvariable $this-colormap4-name
        pack $colormap4.name.label -side left 
        pack $colormap4.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $colormap4.options.replace -variable $this-replace4 -text "Replace field, if it already exists"
        pack $colormap4.options.replace -side left

        sci_frame $colormap5.name
        sci_frame $colormap5.options
        pack $colormap5.name $colormap5.options -side top -fill x -expand yes -padx 5p

        sci_label $colormap5.name.label -text "Name"
        sci_entry $colormap5.name.entry -textvariable $this-colormap5-name
        pack $colormap5.name.label -side left 
        pack $colormap5.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $colormap5.options.replace -variable $this-replace5 -text "Replace field, if it already exists"
        pack $colormap5.options.replace -side left

        sci_frame $colormap6.name
        sci_frame $colormap6.options
        pack $colormap6.name $colormap6.options -side top -fill x -expand yes -padx 5p

        sci_label $colormap6.name.label -text "Name"
        sci_entry $colormap6.name.entry -textvariable $this-colormap6-name
        pack $colormap6.name.label -side left 
        pack $colormap6.name.entry -side left -fill x -expand yes

        sci_checkbutton $colormap6.options.replace -variable $this-replace6 -text "Replace field, if it already exists"
        pack $colormap6.options.replace -side left

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
