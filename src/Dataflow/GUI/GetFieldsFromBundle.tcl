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


itcl::class SCIRun_Bundle_GetFieldsFromBundle {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name GetFieldsFromBundle

        initGlobal $this-field1-listbox ""
        initGlobal $this-field2-listbox ""
        initGlobal $this-field3-listbox ""
        initGlobal $this-field4-listbox ""
        initGlobal $this-field5-listbox ""
        initGlobal $this-field6-listbox ""
        initGlobal $this-field1-entry ""
        initGlobal $this-field2-entry ""
        initGlobal $this-field3-entry ""
        initGlobal $this-field4-entry ""
        initGlobal $this-field5-entry ""
        initGlobal $this-field6-entry ""
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
        global $this-field-selection
        global $this-field1-listbox
        global $this-field2-listbox
        global $this-field3-listbox
        global $this-field4-listbox
        global $this-field5-listbox
        global $this-field6-listbox
        global $this-field1-entry
        global $this-field2-entry
        global $this-field3-entry
        global $this-field4-entry
        global $this-field5-entry
        global $this-field6-entry

        sci_toplevel $w 

        wm minsize $w 100 150
        
        sci_labeledframe $w.frame -labeltext "BUNDLE FIELD OUTPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill both -expand yes

        sci_tabnotebook $childframe.pw -height 160 -width 400 -tabpos n
        $childframe.pw add -label "Field1"
        $childframe.pw add -label "Field2" 
        $childframe.pw add -label "Field3" 
        $childframe.pw add -label "Field4"
        $childframe.pw add -label "Field5" 
        $childframe.pw add -label "Field6" 
        $childframe.pw select 0

        pack $childframe.pw -fill both -expand yes

        set field1 [$childframe.pw childsite 0]
        set field2 [$childframe.pw childsite 1]
        set field3 [$childframe.pw childsite 2]
        set field4 [$childframe.pw childsite 3]
        set field5 [$childframe.pw childsite 4]
        set field6 [$childframe.pw childsite 5]

        sci_frame $field1.name
        sci_frame $field1.sel
        pack $field1.name -side top -fill x -expand no -padx 5p
        pack $field1.sel -side top -fill both -expand yes -padx 5p

        sci_label $field1.name.label -text "Name"
        sci_entry $field1.name.entry -textvariable $this-field1-name
        set $this-field1-entry $field1.name.entry
        pack $field1.name.label -side left 
        pack $field1.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $field1.sel.listbox  -selectioncommand [format "%s ChooseField1" $this]
        set $this-field1-listbox $field1.sel.listbox
        $field1.sel.listbox component listbox configure -listvariable $this-field-selection -selectmode browse
        pack $field1.sel.listbox -fill both -expand yes

        sci_frame $field2.name
        sci_frame $field2.sel 
        pack $field2.name -side top -fill x -expand no -padx 5p
        pack $field2.sel -side top -fill both -expand yes -padx 5p

        sci_label $field2.name.label -text "Name"
        sci_entry $field2.name.entry -textvariable $this-field2-name
        set $this-field2-entry $field2.name.entry
        pack $field2.name.label -side left 
        pack $field2.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $field2.sel.listbox  -selectioncommand [format "%s ChooseField2" $this]
        set $this-field2-listbox $field2.sel.listbox
        $field2.sel.listbox component listbox configure -listvariable $this-field-selection -selectmode browse
        pack $field2.sel.listbox -fill both -expand yes
        
        sci_frame $field3.name
        sci_frame $field3.sel
        pack $field3.name -side top -fill x -expand no -padx 5p
        pack $field3.sel -side top -fill both -expand yes -padx 5p

        sci_label $field3.name.label -text "Name"
        sci_entry $field3.name.entry -textvariable $this-field3-name
        set $this-field3-entry $field3.name.entry
        pack $field3.name.label -side left 
        pack $field3.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $field3.sel.listbox  -selectioncommand [format "%s ChooseField3" $this]
        set $this-field3-listbox $field3.sel.listbox
        $field3.sel.listbox component listbox configure -listvariable $this-field-selection -selectmode browse
        pack $field3.sel.listbox -fill both -expand yes

        sci_frame $field4.name
        sci_frame $field4.sel
        pack $field4.name -side top -fill x -expand no -padx 5p
        pack $field4.sel -side top -fill both -expand yes -padx 5p

        sci_label $field4.name.label -text "Name"
        sci_entry $field4.name.entry -textvariable $this-field4-name
        set $this-field4-entry $field4.name.entry
        pack $field4.name.label -side left 
        pack $field4.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $field4.sel.listbox  -selectioncommand [format "%s ChooseField4" $this]
        set $this-field4-listbox $field4.sel.listbox
        $field4.sel.listbox component listbox configure -listvariable $this-field-selection -selectmode browse
        pack $field4.sel.listbox -fill both -expand yes

        sci_frame $field5.name
        sci_frame $field5.sel
        pack $field5.name -side top -fill x -expand no -padx 5p
        pack $field5.sel -side top -fill both -expand yes -padx 5p

        label $field5.name.label -text "Name"
        entry $field5.name.entry -textvariable $this-field5-name
        set $this-field5-entry $field5.name.entry
        pack $field5.name.label -side left 
        pack $field5.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $field5.sel.listbox  -selectioncommand [format "%s ChooseField5" $this]
        set $this-field5-listbox $field5.sel.listbox
        $field5.sel.listbox component listbox configure -listvariable $this-field-selection -selectmode browse
        pack $field5.sel.listbox -fill both -expand yes
        
        sci_frame $field6.name
        sci_frame $field6.sel
        pack $field6.name -side top -fill x -expand no -padx 5p
        pack $field6.sel -side top -fill both -expand yes -padx 5p

        sci_label $field6.name.label -text "Name"
        sci_entry $field6.name.entry -textvariable $this-field6-name
        set $this-field6-entry $field6.name.entry
        pack $field6.name.label -side left 
        pack $field6.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $field6.sel.listbox  -selectioncommand [format "%s ChooseField6" $this]
        set $this-field6-listbox $field6.sel.listbox
        $field6.sel.listbox component listbox configure -listvariable $this-field-selection -selectmode browse
        pack $field6.sel.listbox -fill both -expand yes

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
    
    
    method ChooseField1 { } {
        global $this-field1-listbox
        global $this-field1-name
        global $this-field-selection
        
        set fieldnum [[set $this-field1-listbox] curselection]
        if [expr [string equal $fieldnum ""] == 0] {
            set $this-field1-name  [lindex [set $this-field-selection] $fieldnum] 
        }
    }

    method ChooseField2 { } {
        global $this-field2-listbox
        global $this-field2-name
        global $this-field-selection
        
        set fieldnum [[set $this-field2-listbox] curselection]
        if [expr [string equal $fieldnum ""] == 0] {
            set $this-field2-name  [lindex [set $this-field-selection] $fieldnum] 
        }
    }

    method ChooseField3 { } {
        global $this-field3-listbox
        global $this-field3-name
        global $this-field-selection
        
        set fieldnum [[set $this-field3-listbox] curselection]
        if [expr [string equal $fieldnum ""] == 0] {
            set $this-field3-name  [lindex [set $this-field-selection] $fieldnum] 
        }
    }

    method ChooseField4 { } {
        global $this-field4-listbox
        global $this-field4-name
        global $this-field-selection
        
        set fieldnum [[set $this-field4-listbox] curselection]
        if [expr [string equal $fieldnum ""] == 0] {
            set $this-field4-name  [lindex [set $this-field-selection] $fieldnum] 
        }
    }

    method ChooseField5 { } {
        global $this-field5-listbox
        global $this-field5-name
        global $this-field-selection
        
        set fieldnum [[set $this-field5-listbox] curselection]
        if [expr [string equal $fieldnum ""] == 0] {
            set $this-field5-name  [lindex [set $this-field-selection] $fieldnum] 
        }
    }

    method ChooseField6 { } {
        global $this-field6-listbox
        global $this-field6-name
        global $this-field-selection
        
        set fieldnum [[set $this-field6-listbox] curselection]
        if [expr [string equal $fieldnum ""] == 0] {
            set $this-field6-name  [lindex [set $this-field-selection] $fieldnum] 
        }
    }
    
}
