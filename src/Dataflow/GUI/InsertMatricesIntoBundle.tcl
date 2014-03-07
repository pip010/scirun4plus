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


itcl::class SCIRun_Bundle_InsertMatricesIntoBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name InsertMatricesIntoBundle
    }

    method ui {} {
        
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        # input matrix names

        global $this-matrix1-name
        global $this-matrix2-name
        global $this-matrix3-name
        global $this-matrix4-name
        global $this-matrix5-name
        global $this-matrix6-name
        global $this-bundlename

        toplevel $w 

        wm minsize $w 100 150

        
        sci_labeledframe $w.frame -labeltext "MATRIX INPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill x
        sci_frame $w.frame2
        pack $w.frame2 -fill x
        sci_label $w.frame2.label -text "Name bundle object :"
        sci_entry $w.frame2.entry -textvariable $this-bundlename
        pack $w.frame2.label -side left 
        pack $w.frame2.entry -side left -fill x

        sci_tabnotebook $childframe.pw -height 100 -width 450 -tabpos n
        $childframe.pw add -label "Matrix1"
        $childframe.pw add -label "Matrix2" 
        $childframe.pw add -label "Matrix3" 
        $childframe.pw add -label "Matrix4"
        $childframe.pw add -label "Matrix5" 
        $childframe.pw add -label "Matrix6" 
        $childframe.pw select 0

        pack $childframe.pw -fill x -expand yes

        set matrix1 [$childframe.pw childsite 0]
        set matrix2 [$childframe.pw childsite 1]
        set matrix3 [$childframe.pw childsite 2]
        set matrix4 [$childframe.pw childsite 3]
        set matrix5 [$childframe.pw childsite 4]
        set matrix6 [$childframe.pw childsite 5]

        sci_frame $matrix1.name
        sci_frame $matrix1.options
        pack $matrix1.name $matrix1.options -side top -fill x -expand yes -padx 5p

        sci_label $matrix1.name.label -text "Name"
        sci_entry $matrix1.name.entry -textvariable $this-matrix1-name
        pack $matrix1.name.label -side left 
        pack $matrix1.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $matrix1.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $matrix1.options.replace -side left

        sci_frame $matrix2.name
        sci_frame $matrix2.options
        pack $matrix2.name $matrix2.options -side top -fill x -expand yes -padx 5p

        sci_label $matrix2.name.label -text "Name"
        sci_entry $matrix2.name.entry -textvariable $this-matrix2-name
        pack $matrix2.name.label -side left 
        pack $matrix2.name.entry -side left -fill x -expand yes

        sci_checkbutton $matrix2.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $matrix2.options.replace -side left
        
        sci_frame $matrix3.name
        sci_frame $matrix3.options
        pack $matrix3.name $matrix3.options -side top -fill x -expand yes -padx 5p

        sci_label $matrix3.name.label -text "Name"
        sci_entry $matrix3.name.entry -textvariable $this-matrix3-name
        pack $matrix3.name.label -side left 
        pack $matrix3.name.entry -side left -fill x -expand yes

        checkbutton $matrix3.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $matrix3.options.replace -side left

        sci_frame $matrix4.name
        sci_frame $matrix4.options
        pack $matrix4.name $matrix4.options -side top -fill x -expand yes -padx 5p

        sci_label $matrix4.name.label -text "Name"
        sci_entry $matrix4.name.entry -textvariable $this-matrix4-name
        pack $matrix4.name.label -side left 
        pack $matrix4.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $matrix4.options.replace -variable $this-replace4 -text "Replace field, if it already exists"
        pack $matrix4.options.replace -side left

        sci_frame $matrix5.name
        sci_frame $matrix5.options
        pack $matrix5.name $matrix5.options -side top -fill x -expand yes -padx 5p

        sci_label $matrix5.name.label -text "Name"
        sci_entry $matrix5.name.entry -textvariable $this-matrix5-name
        pack $matrix5.name.label -side left 
        pack $matrix5.name.entry -side left -fill x -expand yes

        sci_checkbutton $matrix5.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $matrix5.options.replace -side left
        
        sci_frame $matrix6.name
        sci_frame $matrix6.options
        pack $matrix6.name $matrix6.options -side top -fill x -expand yes -padx 5p

        sci_label $matrix6.name.label -text "Name"
        sci_entry $matrix6.name.entry -textvariable $this-matrix6-name
        pack $matrix6.name.label -side left 
        pack $matrix6.name.entry -side left -fill x -expand yes

        sci_checkbutton $matrix6.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $matrix6.options.replace -side left
        
        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
