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


itcl::class SCIRun_Bundle_GetStringsFromBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GetStringsFromBundle

        initGlobal $this-string1-listbox ""
        initGlobal $this-string2-listbox ""
        initGlobal $this-string3-listbox ""
        initGlobal $this-string4-listbox ""
        initGlobal $this-string5-listbox ""
        initGlobal $this-string6-listbox ""
        initGlobal $this-string1-entry ""
        initGlobal $this-string2-entry ""
        initGlobal $this-string3-entry ""
        initGlobal $this-string4-entry ""
        initGlobal $this-string5-entry ""
        initGlobal $this-string6-entry ""
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
        global $this-string-selection
        global $this-string1-listbox
        global $this-string2-listbox
        global $this-string3-listbox
        global $this-string4-listbox
        global $this-string5-listbox
        global $this-string6-listbox
        global $this-string1-entry
        global $this-string2-entry
        global $this-string3-entry
        global $this-string4-entry
        global $this-string5-entry
        global $this-string6-entry

        sci_toplevel $w 

        wm minsize $w 100 150
        
        sci_labeledframe $w.frame -labeltext "STRING OUTPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill both -expand yes

        sci_tabnotebook $childframe.pw -width 400 -tabpos n
        $childframe.pw add -label "String1"
        $childframe.pw add -label "String2" 
        $childframe.pw add -label "String3" 
        $childframe.pw add -label "String4"
        $childframe.pw add -label "String5" 
        $childframe.pw add -label "String6" 
        $childframe.pw select 0

        pack $childframe.pw -fill both -expand yes

        set string1 [$childframe.pw childsite 0]
        set string2 [$childframe.pw childsite 1]
        set string3 [$childframe.pw childsite 2]
        set string4 [$childframe.pw childsite 3]
        set string5 [$childframe.pw childsite 4]
        set string6 [$childframe.pw childsite 5]

        sci_frame $string1.name
        sci_frame $string1.sel
        pack $string1.name -side top -fill x -expand no -padx 5p
        pack $string1.sel -side top -fill both -expand yes -padx 5p

        sci_label $string1.name.label -text "Name"
        sci_entry $string1.name.entry -textvariable $this-string1-name
        set $this-string1-entry $string1.name.entry
        pack $string1.name.label -side left 
        pack $string1.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $string1.sel.listbox  -selectioncommand [format "%s ChooseString1" $this]
        set $this-string1-listbox $string1.sel.listbox
        $string1.sel.listbox component listbox configure -listvariable $this-string-selection -selectmode browse
        pack $string1.sel.listbox -fill both -expand yes

        sci_frame $string2.name
        sci_frame $string2.sel
        pack $string2.name -side top -fill x -expand no -padx 5p
        pack $string2.sel -side top -fill both -expand yes -padx 5p

        sci_label $string2.name.label -text "Name"
        sci_entry $string2.name.entry -textvariable $this-string2-name
        set $this-string2-entry $string2.name.entry
        pack $string2.name.label -side left 
        pack $string2.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $string2.sel.listbox  -selectioncommand [format "%s ChooseString2" $this]
        set $this-string2-listbox $string2.sel.listbox
        $string2.sel.listbox component listbox configure -listvariable $this-string-selection -selectmode browse
        pack $string2.sel.listbox -fill both -expand yes
        
        sci_frame $string3.name
        sci_frame $string3.sel
        pack $string3.name -side top -fill x -expand no -padx 5p
        pack $string3.sel -side top -fill both -expand yes -padx 5p

        sci_label $string3.name.label -text "Name"
        sci_entry $string3.name.entry -textvariable $this-string3-name
        set $this-string3-entry $string3.name.entry
        pack $string3.name.label -side left 
        pack $string3.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $string3.sel.listbox  -selectioncommand [format "%s ChooseString3" $this]
        set $this-string3-listbox $string3.sel.listbox
        $string3.sel.listbox component listbox configure -listvariable $this-string-selection -selectmode browse
        pack $string3.sel.listbox -fill both -expand yes

        sci_frame $string4.name
        sci_frame $string4.sel
        pack $string4.name -side top -fill x -expand no -padx 5p
        pack $string4.sel -side top -fill both -expand yes -padx 5p

        sci_label $string4.name.label -text "Name"
        sci_entry $string4.name.entry -textvariable $this-string4-name
        set $this-string4-entry $string4.name.entry
        pack $string4.name.label -side left 
        pack $string4.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $string4.sel.listbox  -selectioncommand [format "%s ChooseString4" $this]
        set $this-string4-listbox $string4.sel.listbox
        $string4.sel.listbox component listbox configure -listvariable $this-string-selection -selectmode browse
        pack $string4.sel.listbox -fill both -expand yes

        sci_frame $string5.name
        sci_frame $string5.sel
        pack $string5.name -side top -fill x -expand no -padx 5p
        pack $string5.sel -side top -fill both -expand yes -padx 5p

        sci_label $string5.name.label -text "Name"
        sci_entry $string5.name.entry -textvariable $this-string5-name
        set $this-string5-entry $string5.name.entry
        pack $string5.name.label -side left 
        pack $string5.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $string5.sel.listbox  -selectioncommand [format "%s ChooseString5" $this]
        set $this-string5-listbox $string5.sel.listbox
        $string5.sel.listbox component listbox configure -listvariable $this-string-selection -selectmode browse
        pack $string5.sel.listbox -fill both -expand yes
        
        sci_frame $string6.name
        sci_frame $string6.sel
        pack $string6.name -side top -fill x -expand no -padx 5p
        pack $string6.sel -side top -fill both -expand yes -padx 5p

        sci_label $string6.name.label -text "Name"
        sci_entry $string6.name.entry -textvariable $this-string6-name
        set $this-string6-entry $string6.name.entry
        pack $string6.name.label -side left 
        pack $string6.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $string6.sel.listbox  -selectioncommand [format "%s ChooseString6" $this]
        set $this-string6-listbox $string6.sel.listbox
        $string6.sel.listbox component listbox configure -listvariable $this-string-selection -selectmode browse
        pack $string6.sel.listbox -fill both -expand yes


        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
    
    
    method ChooseString1 { } {
        global $this-string1-listbox
        global $this-string1-name
        global $this-string-selection
        
        set stringnum [[set $this-string1-listbox] curselection]
        if [expr [string equal $stringnum ""] == 0] {
            set $this-string1-name  [lindex [set $this-string-selection] $stringnum] 
        }
    }

    method ChooseString2 { } {
        global $this-string2-listbox
        global $this-string2-name
        global $this-string-selection
        
        set stringnum [[set $this-string2-listbox] curselection]
        if [expr [string equal $stringnum ""] == 0] {
            set $this-string2-name  [lindex [set $this-string-selection] $stringnum] 
        }
    }

    method ChooseString3 { } {
        global $this-string3-listbox
        global $this-string3-name
        global $this-string-selection
        
        set stringnum [[set $this-string3-listbox] curselection]
        if [expr [string equal $stringnum ""] == 0] {
            set $this-string3-name  [lindex [set $this-string-selection] $stringnum] 
        }
    }

    method ChooseString4 { } {
        global $this-string4-listbox
        global $this-string4-name
        global $this-string-selection
        
        set stringnum [[set $this-string4-listbox] curselection]
        if [expr [string equal $stringnum ""] == 0] {
            set $this-string4-name  [lindex [set $this-string-selection] $stringnum] 
        }
    }

    method ChooseString5 { } {
        global $this-string5-listbox
        global $this-string5-name
        global $this-string-selection
        
        set stringnum [[set $this-string5-listbox] curselection]
        if [expr [string equal $stringnum ""] == 0] {
            set $this-string5-name  [lindex [set $this-string-selection] $stringnum] 
        }
    }

    method ChooseString6 { } {
        global $this-string6-listbox
        global $this-string6-name
        global $this-string-selection
        
        set stringnum [[set $this-string6-listbox] curselection]
        if [expr [string equal $stringnum ""] == 0] {
            set $this-string6-name  [lindex [set $this-string-selection] $stringnum] 
        }
    }
    
}
