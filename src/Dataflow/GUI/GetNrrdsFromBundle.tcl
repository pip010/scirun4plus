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


itcl::class SCIRun_Bundle_GetNrrdsFromBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GetNrrdsFromBundle

        initGlobal $this-nrrd1-listbox ""
        initGlobal $this-nrrd2-listbox ""
        initGlobal $this-nrrd3-listbox ""
        initGlobal $this-nrrd4-listbox ""
        initGlobal $this-nrrd5-listbox ""
        initGlobal $this-nrrd6-listbox ""
        initGlobal $this-nrrd1-entry ""
        initGlobal $this-nrrd2-entry ""
        initGlobal $this-nrrd3-entry ""
        initGlobal $this-nrrd4-entry ""
        initGlobal $this-nrrd5-entry ""
        initGlobal $this-nrrd6-entry ""
    }

    method ui {} {
        
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        # input nrrd names

        global $this-nrrd1-name
        global $this-nrrd2-name
        global $this-nrrd3-name
        global $this-nrrd4-name
        global $this-nrrd5-name
        global $this-nrrd6-name
        global $this-nrrd-selection
        global $this-nrrd1-listbox
        global $this-nrrd2-listbox
        global $this-nrrd3-listbox
        global $this-nrrd4-listbox
        global $this-nrrd5-listbox
        global $this-nrrd6-listbox
        global $this-nrrd1-entry
        global $this-nrrd2-entry
        global $this-nrrd3-entry
        global $this-nrrd4-entry
        global $this-nrrd5-entry
        global $this-nrrd6-entry
        global $this-transposenrrd1
        global $this-transposenrrd2
        global $this-transposenrrd3
        global $this-transposenrrd4
        global $this-transposenrrd5
        global $this-transposenrrd6
        
        sci_toplevel $w 

        wm minsize $w 100 150

        
        sci_labeledframe $w.frame -labeltext "BUNDLE NRRD OUTPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill both -expand yes

        sci_tabnotebook $childframe.pw -height 160 -width 400 -tabpos n
        $childframe.pw add -label "Nrrd1"
        $childframe.pw add -label "Nrrd2" 
        $childframe.pw add -label "Nrrd3" 
        $childframe.pw add -label "Nrrd4"
        $childframe.pw add -label "Nrrd5" 
        $childframe.pw add -label "Nrrd6" 
        $childframe.pw select 0

        pack $childframe.pw -fill both -expand yes

        set nrrd1 [$childframe.pw childsite 0]
        set nrrd2 [$childframe.pw childsite 1]
        set nrrd3 [$childframe.pw childsite 2]
        set nrrd4 [$childframe.pw childsite 3]
        set nrrd5 [$childframe.pw childsite 4]
        set nrrd6 [$childframe.pw childsite 5]

        sci_frame $nrrd1.name
        sci_frame $nrrd1.sel
        sci_frame $nrrd1.transpose
        pack $nrrd1.name -side top -fill x -padx 5p
        pack $nrrd1.transpose -side top -fill x -padx 5p
        pack $nrrd1.sel -side top -fill both -expand yes -padx 5p

        sci_label $nrrd1.name.label -text "Name"
        sci_entry $nrrd1.name.entry -textvariable $this-nrrd1-name
        set $this-nrrd1-entry $nrrd1.name.entry
        pack $nrrd1.name.label -side left 
        pack $nrrd1.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $nrrd1.sel.listbox  -selectioncommand [format "%s ChooseNrrd1" $this] 
        set $this-nrrd1-listbox $nrrd1.sel.listbox
        $nrrd1.sel.listbox component listbox configure -listvariable $this-nrrd-selection -selectmode browse
        pack $nrrd1.sel.listbox -fill both -expand yes
        sci_checkbutton $nrrd1.transpose.cb -variable $this-transposenrrd1 -text "Assume matrix data is transposed"
        pack $nrrd1.transpose.cb -side left -fill x

        sci_frame $nrrd2.name
        sci_frame $nrrd2.sel
        sci_frame $nrrd2.transpose    
        pack $nrrd2.name -side top -fill x -padx 5p
        pack $nrrd2.transpose -side top -fill x -padx 5p
        pack $nrrd2.sel -side top -fill both -expand yes -padx 5p

        sci_label $nrrd2.name.label -text "Name"
        sci_entry $nrrd2.name.entry -textvariable $this-nrrd2-name
        set $this-nrrd2-entry $nrrd2.name.entry
        pack $nrrd2.name.label -side left 
        pack $nrrd2.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $nrrd2.sel.listbox  -selectioncommand [format "%s ChooseNrrd2" $this] 
        set $this-nrrd2-listbox $nrrd2.sel.listbox
        $nrrd2.sel.listbox component listbox configure -listvariable $this-nrrd-selection -selectmode browse
        pack $nrrd2.sel.listbox -fill both -expand yes
        sci_checkbutton $nrrd2.transpose.cb -variable $this-transposenrrd1 -text "Assume matrix data is transposed"
        pack $nrrd2.transpose.cb -side left -fill x
        
        sci_frame $nrrd3.name
        sci_frame $nrrd3.sel
        sci_frame $nrrd3.transpose
        pack $nrrd3.name -side top -fill x -padx 5p
        pack $nrrd3.transpose -side top -fill x -padx 5p
        pack $nrrd3.sel -side top -fill both -expand yes -padx 5p

        sci_label $nrrd3.name.label -text "Name"
        sci_entry $nrrd3.name.entry -textvariable $this-nrrd3-name
        set $this-nrrd3-entry $nrrd3.name.entry
        pack $nrrd3.name.label -side left 
        pack $nrrd3.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $nrrd3.sel.listbox  -selectioncommand [format "%s ChooseNrrd3" $this] 
        set $this-nrrd3-listbox $nrrd3.sel.listbox
        $nrrd3.sel.listbox component listbox configure -listvariable $this-nrrd-selection -selectmode browse
        pack $nrrd3.sel.listbox -fill both -expand yes
        sci_checkbutton $nrrd3.transpose.cb -variable $this-transposenrrd1 -text "Assume matrix data is transposed"
        pack $nrrd3.transpose.cb -side left -fill x

        sci_frame $nrrd4.name
        sci_frame $nrrd4.sel
        sci_frame $nrrd4.transpose
        pack $nrrd4.name -side top -fill x -padx 5p
        pack $nrrd4.transpose -side top -fill x -padx 5p
        pack $nrrd4.sel -side top -fill both -expand yes -padx 5p

        sci_label $nrrd4.name.label -text "Name"
        sci_entry $nrrd4.name.entry -textvariable $this-nrrd4-name
        set $this-nrrd4-entry $nrrd4.name.entry
        pack $nrrd4.name.label -side left 
        pack $nrrd4.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $nrrd4.sel.listbox  -selectioncommand [format "%s ChooseNrrd4" $this] 
        set $this-nrrd4-listbox $nrrd4.sel.listbox
        $nrrd4.sel.listbox component listbox configure -listvariable $this-nrrd-selection -selectmode browse
        pack $nrrd4.sel.listbox -fill both -expand yes
        sci_checkbutton $nrrd4.transpose.cb -variable $this-transposenrrd4 -text "Assume matrix data is transposed"
        pack $nrrd4.transpose.cb -side left -fill x

        sci_frame $nrrd5.name
        sci_frame $nrrd5.sel
        sci_frame $nrrd5.transpose    
        pack $nrrd5.name -side top -fill x -padx 5p
        pack $nrrd5.transpose -side top -fill x -padx 5p
        pack $nrrd5.sel -side top -fill both -expand yes -padx 5p

        sci_label $nrrd5.name.label -text "Name"
        sci_entry $nrrd5.name.entry -textvariable $this-nrrd5-name
        set $this-nrrd5-entry $nrrd5.name.entry
        pack $nrrd5.name.label -side left 
        pack $nrrd5.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $nrrd5.sel.listbox  -selectioncommand [format "%s ChooseNrrd5" $this] 
        set $this-nrrd5-listbox $nrrd5.sel.listbox
        $nrrd5.sel.listbox component listbox configure -listvariable $this-nrrd-selection -selectmode browse
        pack $nrrd5.sel.listbox -fill both -expand yes
        sci_checkbutton $nrrd5.transpose.cb -variable $this-transposenrrd1 -text "Assume matrix data is transposed"
        pack $nrrd5.transpose.cb -side left -fill x
        
        sci_frame $nrrd6.name
        sci_frame $nrrd6.sel
        sci_frame $nrrd6.transpose
        pack $nrrd6.name -side top -fill x -padx 5p
        pack $nrrd6.transpose -side top -fill x -padx 5p
        pack $nrrd6.sel -side top -fill both -expand yes -padx 5p

        sci_label $nrrd6.name.label -text "Name"
        sci_entry $nrrd6.name.entry -textvariable $this-nrrd6-name
        set $this-nrrd6-entry $nrrd6.name.entry
        pack $nrrd6.name.label -side left 
        pack $nrrd6.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $nrrd6.sel.listbox  -selectioncommand [format "%s ChooseNrrd6" $this] 
        set $this-nrrd6-listbox $nrrd6.sel.listbox
        $nrrd6.sel.listbox component listbox configure -listvariable $this-nrrd-selection -selectmode browse
        pack $nrrd6.sel.listbox -fill both -expand yes
        sci_checkbutton $nrrd6.transpose.cb -variable $this-transposenrrd1 -text "Assume matrix data is transposed"
        pack $nrrd6.transpose.cb -side left -fill x

        makeSciButtonPanel $w $w $this
    }
    
    
    method ChooseNrrd1 { } {
        global $this-nrrd1-listbox
        global $this-nrrd1-name
        global $this-nrrd-selection
        
        set nrrdnum [[set $this-nrrd1-listbox] curselection]
        if [expr [string equal $nrrdnum ""] == 0] {
            set $this-nrrd1-name  [lindex [set $this-nrrd-selection] $nrrdnum] 
        }
    }

    method ChooseNrrd2 { } {
        global $this-nrrd2-listbox
        global $this-nrrd2-name
        global $this-nrrd-selection
        
        set nrrdnum [[set $this-nrrd2-listbox] curselection]
        if [expr [string equal $nrrdnum ""] == 0] {
            set $this-nrrd2-name  [lindex [set $this-nrrd-selection] $nrrdnum] 
        }
    }

    method ChooseNrrd3 { } {
        global $this-nrrd3-listbox
        global $this-nrrd3-name
        global $this-nrrd-selection
        
        set nrrdnum [[set $this-nrrd3-listbox] curselection]
        if [expr [string equal $nrrdnum ""] == 0] {
            set $this-nrrd3-name  [lindex [set $this-nrrd-selection] $nrrdnum] 
        }
    }

    method ChooseNrrd4 { } {
        global $this-nrrd4-listbox
        global $this-nrrd4-name
        global $this-nrrd-selection
        
        set nrrdnum [[set $this-nrrd4-listbox] curselection]
        if [expr [string equal $nrrdnum ""] == 0] {
            set $this-nrrd4-name  [lindex [set $this-nrrd-selection] $nrrdnum] 
        }
    }

    method ChooseNrrd5 { } {
        global $this-nrrd5-listbox
        global $this-nrrd5-name
        global $this-nrrd-selection
        
        set nrrdnum [[set $this-nrrd5-listbox] curselection]
        if [expr [string equal $nrrdnum ""] == 0] {
            set $this-nrrd5-name  [lindex [set $this-nrrd-selection] $nrrdnum] 
        }
    }

    method ChooseNrrd6 { } {
        global $this-nrrd6-listbox
        global $this-nrrd6-name
        global $this-nrrd-selection
        
        set nrrdnum [[set $this-nrrd6-listbox] curselection]
        if [expr [string equal $nrrdnum ""] == 0] {
            set $this-nrrd6-name  [lindex [set $this-nrrd-selection] $nrrdnum] 
        }
    }
    
}
