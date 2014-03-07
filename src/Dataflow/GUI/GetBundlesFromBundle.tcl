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


itcl::class SCIRun_Bundle_GetBundlesFromBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GetBundlesFromBundle

        initGlobal $this-bundle1-listbox ""
        initGlobal $this-bundle2-listbox ""
        initGlobal $this-bundle3-listbox ""
        initGlobal $this-bundle4-listbox ""
        initGlobal $this-bundle5-listbox ""
        initGlobal $this-bundle6-listbox ""
        initGlobal $this-bundle1-entry ""
        initGlobal $this-bundle2-entry ""
        initGlobal $this-bundle3-entry ""
        initGlobal $this-bundle4-entry ""
        initGlobal $this-bundle5-entry ""
        initGlobal $this-bundle6-entry ""
    }

    method ui {} {
        
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        # input bundle names

        global $this-bundle1-name
        global $this-bundle2-name
        global $this-bundle3-name
        global $this-bundle4-name
        global $this-bundle5-name
        global $this-bundle6-name
        global $this-bundle-selection
        global $this-bundle1-listbox
        global $this-bundle2-listbox
        global $this-bundle3-listbox
        global $this-bundle4-listbox
        global $this-bundle5-listbox
        global $this-bundle6-listbox
        global $this-bundle1-entry
        global $this-bundle2-entry
        global $this-bundle3-entry
        global $this-bundle4-entry
        global $this-bundle5-entry
        global $this-bundle6-entry

        sci_toplevel $w 

        wm minsize $w 100 150

        
        sci_labeledframe $w.frame -labeltext "BUNDLE SUB-BUNDLE OUTPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill both -expand yes

        sci_tabnotebook $childframe.pw -width 400 -tabpos n
        $childframe.pw add -label "Bundle1"
        $childframe.pw add -label "Bundle2" 
        $childframe.pw add -label "Bundle3" 
        $childframe.pw add -label "Bundle4"
        $childframe.pw add -label "Bundle5" 
        $childframe.pw add -label "Bundle6" 
        $childframe.pw select 0

        pack $childframe.pw -fill both -expand yes

        set bundle1 [$childframe.pw childsite 0]
        set bundle2 [$childframe.pw childsite 1]
        set bundle3 [$childframe.pw childsite 2]
        set bundle4 [$childframe.pw childsite 3]
        set bundle5 [$childframe.pw childsite 4]
        set bundle6 [$childframe.pw childsite 5]

        sci_frame $bundle1.name
        sci_frame $bundle1.sel
        pack $bundle1.name -side top -fill x -expand no -padx 5p
        pack $bundle1.sel -side top -fill both -expand yes -padx 5p

        sci_label $bundle1.name.label -text "Name"
        sci_entry $bundle1.name.entry -textvariable $this-bundle1-name
        set $this-bundle1-entry $bundle1.name.entry
        pack $bundle1.name.label -side left
        pack $bundle1.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $bundle1.sel.listbox  -selectioncommand [format "%s ChooseBundle1" $this]
        set $this-bundle1-listbox $bundle1.sel.listbox
        $bundle1.sel.listbox component listbox configure -listvariable $this-bundle-selection -selectmode browse
        pack $bundle1.sel.listbox -fill both -expand yes

        sci_frame $bundle2.name
        sci_frame $bundle2.sel
        pack $bundle2.name -side top -fill x -expand no -padx 5p
        pack $bundle2.sel -side top -fill both -expand yes -padx 5p

        sci_label $bundle2.name.label -text "Name"
        sci_entry $bundle2.name.entry -textvariable $this-bundle2-name
        set $this-bundle2-entry $bundle2.name.entry
        pack $bundle2.name.label -side left 
        pack $bundle2.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $bundle2.sel.listbox  -selectioncommand [format "%s ChooseBundle2" $this]
        set $this-bundle2-listbox $bundle2.sel.listbox
        $bundle2.sel.listbox component listbox configure -listvariable $this-bundle-selection -selectmode browse
        pack $bundle2.sel.listbox -fill both -expand yes
        
        sci_frame $bundle3.name
        sci_frame $bundle3.sel
        pack $bundle3.name -side top -fill x -expand no -padx 5p
        pack $bundle3.sel -side top -fill both -expand yes -padx 5p

        sci_label $bundle3.name.label -text "Name"
        sci_entry $bundle3.name.entry -textvariable $this-bundle3-name
        set $this-bundle3-entry $bundle3.name.entry
        pack $bundle3.name.label -side left 
        pack $bundle3.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $bundle3.sel.listbox  -selectioncommand [format "%s ChooseBundle3" $this]
        set $this-bundle3-listbox $bundle3.sel.listbox
        $bundle3.sel.listbox component listbox configure -listvariable $this-bundle-selection -selectmode browse
        pack $bundle3.sel.listbox -fill both -expand yes

        sci_frame $bundle4.name
        sci_frame $bundle4.sel
        pack $bundle4.name -side top -fill x -expand no -padx 5p
        pack $bundle4.sel -side top -fill both -expand yes -padx 5p

        sci_label $bundle4.name.label -text "Name"
        sci_entry $bundle4.name.entry -textvariable $this-bundle4-name
        set $this-bundle4-entry $bundle4.name.entry
        pack $bundle4.name.label -side left 
        pack $bundle4.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $bundle4.sel.listbox  -selectioncommand [format "%s ChooseBundle4" $this]
        set $this-bundle4-listbox $bundle4.sel.listbox
        $bundle4.sel.listbox component listbox configure -listvariable $this-bundle-selection -selectmode browse
        pack $bundle4.sel.listbox -fill both -expand yes

        sci_frame $bundle5.name
        sci_frame $bundle5.sel
        pack $bundle5.name -side top -fill x -expand no -padx 5p
        pack $bundle5.sel -side top -fill both -expand yes -padx 5p

        sci_label $bundle5.name.label -text "Name"
        sci_entry $bundle5.name.entry -textvariable $this-bundle5-name
        set $this-bundle5-entry $bundle5.name.entry
        pack $bundle5.name.label -side left 
        pack $bundle5.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $bundle5.sel.listbox  -selectioncommand [format "%s ChooseBundle5" $this]
        set $this-bundle5-listbox $bundle5.sel.listbox
        $bundle5.sel.listbox component listbox configure -listvariable $this-bundle-selection -selectmode browse
        pack $bundle5.sel.listbox -fill both -expand yes
        
        sci_frame $bundle6.name
        sci_frame $bundle6.sel
        pack $bundle6.name -side top -fill x -expand no -padx 5p
        pack $bundle6.sel -side top -fill both -expand yes -padx 5p

        sci_label $bundle6.name.label -text "Name"
        sci_entry $bundle6.name.entry -textvariable $this-bundle6-name
        set $this-bundle6-entry $bundle6.name.entry
        pack $bundle6.name.label -side left 
        pack $bundle6.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $bundle6.sel.listbox  -selectioncommand [format "%s ChooseBundle6" $this]
        set $this-bundle6-listbox $bundle6.sel.listbox
        $bundle6.sel.listbox component listbox configure -listvariable $this-bundle-selection -selectmode browse
        pack $bundle6.sel.listbox -fill both -expand yes
        
        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
    
    
    method ChooseBundle1 { } {
        global $this-bundle1-listbox
        global $this-bundle1-name
        global $this-bundle-selection
        
        set bundlenum [[set $this-bundle1-listbox] curselection]
        if [expr [string equal $bundlenum ""] == 0] {
            set $this-bundle1-name  [lindex [set $this-bundle-selection] $bundlenum] 
        }
    }

    method ChooseBundle2 { } {
        global $this-bundle2-listbox
        global $this-bundle2-name
        global $this-bundle-selection
        
        set bundlenum [[set $this-bundle2-listbox] curselection]
        if [expr [string equal $bundlenum ""] == 0] {
            set $this-bundle2-name  [lindex [set $this-bundle-selection] $bundlenum] 
        }
    }

    method ChooseBundle3 { } {
        global $this-bundle3-listbox
        global $this-bundle3-name
        global $this-bundle-selection
        
        set bundlenum [[set $this-bundle3-listbox] curselection]
        if [expr [string equal $bundlenum ""] == 0] {
            set $this-bundle3-name  [lindex [set $this-bundle-selection] $bundlenum] 
        }
    }

    method ChooseBundle4 { } {
        global $this-bundle4-listbox
        global $this-bundle4-name
        global $this-bundle-selection
        
        set bundlenum [[set $this-bundle4-listbox] curselection]
        if [expr [string equal $bundlenum ""] == 0] {
            set $this-bundle4-name  [lindex [set $this-bundle-selection] $bundlenum] 
        }
    }

    method ChooseBundle5 { } {
        global $this-bundle5-listbox
        global $this-bundle5-name
        global $this-bundle-selection
        
        set bundlenum [[set $this-bundle5-listbox] curselection]
        if [expr [string equal $bundlenum ""] == 0] {
            set $this-bundle5-name  [lindex [set $this-bundle-selection] $bundlenum] 
        }
    }

    method ChooseBundle6 { } {
        global $this-bundle6-listbox
        global $this-bundle6-name
        global $this-bundle-selection
        
        set bundlenum [[set $this-bundle6-listbox] curselection]
        if [expr [string equal $bundlenum ""] == 0] {
            set $this-bundle6-name  [lindex [set $this-bundle-selection] $bundlenum] 
        }
    }    
}
