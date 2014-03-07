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


itcl::class SCIRun_Bundle_InsertStringsIntoBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name InsertStringsIntoBundle
    }

    method ui {} {
        
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        # input matrix names

        global $this-string1-name
        global $this-string2-name
        global $this-string3-name
        global $this-string4-name
        global $this-string5-name
        global $this-string6-name
        global $this-bundlename

        sci_toplevel $w 

        wm minsize $w 100 150

        
        sci_labeledframe $w.frame -labeltext "STRING INPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill x
        sci_frame $w.frame2
        pack $w.frame2 -fill x
        sci_label $w.frame2.label -text "Name bundle object :"
        sci_entry $w.frame2.entry -textvariable $this-bundlename
        pack $w.frame2.label -side left 
        pack $w.frame2.entry -side left -fill x

        sci_tabnotebook $childframe.pw -height 100 -width 450 -tabpos n
        $childframe.pw add -label "String1"
        $childframe.pw add -label "String2" 
        $childframe.pw add -label "String3" 
        $childframe.pw add -label "String4"
        $childframe.pw add -label "String5" 
        $childframe.pw add -label "String6" 
        $childframe.pw select 0

        pack $childframe.pw -fill x -expand yes

        set string1 [$childframe.pw childsite 0]
        set string2 [$childframe.pw childsite 1]
        set string3 [$childframe.pw childsite 2]
        set string4 [$childframe.pw childsite 3]
        set string5 [$childframe.pw childsite 4]
        set string6 [$childframe.pw childsite 5]

        sci_frame $string1.name
        sci_frame $string1.options
        pack $string1.name $string1.options -side top -fill x -expand yes -padx 5p

        sci_label $string1.name.label -text "Name"
        sci_entry $string1.name.entry -textvariable $this-string1-name
        pack $string1.name.label -side left 
        pack $string1.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $string1.options.replace -variable $this-replace1 -text "Replace field, if it already exists"
        pack $string1.options.replace -side left
        
        sci_frame $string2.name
        sci_frame $string2.options
        pack $string2.name $string2.options -side top -fill x -expand yes -padx 5p

        sci_label $string2.name.label -text "Name"
        sci_entry $string2.name.entry -textvariable $this-string2-name
        pack $string2.name.label -side left 
        pack $string2.name.entry -side left -fill x -expand yes

        sci_checkbutton $string2.options.replace -variable $this-replace2 -text "Replace field, if it already exists"
        pack $string2.options.replace -side left
        
        sci_frame $string3.name
        sci_frame $string3.options
        pack $string3.name $string3.options -side top -fill x -expand yes -padx 5p

        sci_label $string3.name.label -text "Name"
        sci_entry $string3.name.entry -textvariable $this-string3-name
        pack $string3.name.label -side left 
        pack $string3.name.entry -side left -fill x -expand yes

        sci_checkbutton $string3.options.replace -variable $this-replace3 -text "Replace field, if it already exists"
        pack $string3.options.replace -side left

        sci_frame $string4.name
        sci_frame $string4.options
        pack $string4.name $string4.options -side top -fill x -expand yes -padx 5p

        sci_label $string4.name.label -text "Name"
        sci_entry $string4.name.entry -textvariable $this-string4-name
        pack $string4.name.label -side left 
        pack $string4.name.entry -side left -fill x -expand yes
        
        sci_checkbutton $string4.options.replace -variable $this-replace4 -text "Replace field, if it already exists"
        pack $string4.options.replace -side left
        
        sci_frame $string5.name
        sci_frame $string5.options
        pack $string5.name $string5.options -side top -fill x -expand yes -padx 5p

        sci_label $string5.name.label -text "Name"
        sci_entry $string5.name.entry -textvariable $this-string5-name
        pack $string5.name.label -side left 
        pack $string5.name.entry -side left -fill x -expand yes

        sci_checkbutton $string5.options.replace -variable $this-replace5 -text "Replace field, if it already exists"
        pack $string5.options.replace -side left
        
        sci_frame $string6.name
        sci_frame $string6.options
        pack $string6.name $string6.options -side top -fill x -expand yes -padx 5p

        sci_label $string6.name.label -text "Name"
        sci_entry $string6.name.entry -textvariable $this-string6-name
        pack $string6.name.label -side left 
        pack $string6.name.entry -side left -fill x -expand yes

        sci_checkbutton $string6.options.replace -variable $this-replace6 -text "Replace field, if it already exists"
        pack $string6.options.replace -side left

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}
