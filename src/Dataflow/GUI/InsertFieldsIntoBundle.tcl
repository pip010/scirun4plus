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


itcl::class SCIRun_Bundle_InsertFieldsIntoBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name InsertFieldsIntoBundle
    }

    method ui {} {
        
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        # input matrix names

        global $this-field1-name
        global $this-field2-name
        global $this-field3-name
        global $this-field4-name
        global $this-field5-name
        global $this-field6-name
        global $this-bundlename

        sci_toplevel $w 

        wm minsize $w 100 150
        
        sci_labeledframe $w.frame -labeltext "FIELD INPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill x
        sci_frame $w.frame2
        pack $w.frame2 -fill x
        sci_label $w.frame2.label -text "Name bundle object :"
        sci_entry $w.frame2.entry -textvariable $this-bundlename
        pack $w.frame2.label -side left 
        pack $w.frame2.entry -side left -fill x

        sci_tabnotebook $childframe.pw -height 100 -width 450 -tabpos n
        $childframe.pw add -label "Field1"
        $childframe.pw add -label "Field2" 
        $childframe.pw add -label "Field3" 
        $childframe.pw add -label "Field4"
        $childframe.pw add -label "Field5" 
        $childframe.pw add -label "Field6" 
        $childframe.pw select 0

        pack $childframe.pw -fill x -expand yes

        set field1 [$childframe.pw childsite 0]
        set field2 [$childframe.pw childsite 1]
        set field3 [$childframe.pw childsite 2]
        set field4 [$childframe.pw childsite 3]
        set field5 [$childframe.pw childsite 4]
        set field6 [$childframe.pw childsite 5]

        sci_frame $field1.name
        sci_frame $field1.options
        pack $field1.name $field1.options -side top -fill x -expand yes -padx 5p

        sci_label $field1.name.label -text "Name"
        sci_entry $field1.name.entry -textvariable $this-field1-name
        pack $field1.name.label -side left 
        pack $field1.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $field1.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $field1.options.replace -side left
        
        sci_frame $field2.name
        sci_frame $field2.options
        pack $field2.name $field2.options -side top -fill x -expand yes -padx 5p

        sci_label $field2.name.label -text "Name"
        sci_entry $field2.name.entry -textvariable $this-field2-name
        pack $field2.name.label -side left 
        pack $field2.name.entry -side left -fill x -expand yes

        sci_checkbutton $field2.options.replace -variable $this-replace2 -text "Replace field, if it already exists"
        pack $field2.options.replace -side left
        
        sci_frame $field3.name
        sci_frame $field3.options
        pack $field3.name $field3.options -side top -fill x -expand yes -padx 5p

        sci_label $field3.name.label -text "Name"
        sci_entry $field3.name.entry -textvariable $this-field3-name
        pack $field3.name.label -side left 
        pack $field3.name.entry -side left -fill x -expand yes

        sci_checkbutton $field3.options.replace -variable $this-replace3 -text "Replace field, if it already exists"
        pack $field3.options.replace -side left

        sci_frame $field4.name
        sci_frame $field4.options
        pack $field4.name $field4.options -side top -fill x -expand yes -padx 5p

        sci_label $field4.name.label -text "Name"
        sci_entry $field4.name.entry -textvariable $this-field4-name
        pack $field4.name.label -side left 
        pack $field4.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $field4.options.replace -variable $this-replace4 -text "Replace field, if it already exists"
        pack $field4.options.replace -side left
        
        sci_frame $field5.name
        sci_frame $field5.options
        pack $field5.name $field5.options -side top -fill x -expand yes -padx 5p

        sci_label $field5.name.label -text "Name"
        sci_entry $field5.name.entry -textvariable $this-field5-name
        pack $field5.name.label -side left 
        pack $field5.name.entry -side left -fill x -expand yes

        sci_checkbutton $field5.options.replace -variable $this-replace5 -text "Replace field, if it already exists"
        pack $field5.options.replace -side left
        
        sci_frame $field6.name
        sci_frame $field6.options
        pack $field6.name $field6.options -side top -fill x -expand yes -padx 5p

        sci_label $field6.name.label -text "Name"
        sci_entry $field6.name.entry -textvariable $this-field6-name
        pack $field6.name.label -side left 
        pack $field6.name.entry -side left -fill x -expand yes

        sci_checkbutton $field6.options.replace -variable $this-replace6 -text "Replace field, if it already exists"
        pack $field6.options.replace -side left
        
        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
