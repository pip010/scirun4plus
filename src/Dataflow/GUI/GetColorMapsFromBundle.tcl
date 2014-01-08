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


itcl::class SCIRun_Bundle_GetColorMapsFromBundle {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GetColorMapsFromBundle

        initGlobal $this-colormap1-listbox ""
        initGlobal $this-colormap2-listbox ""
        initGlobal $this-colormap3-listbox ""
        initGlobal $this-colormap4-listbox ""
        initGlobal $this-colormap5-listbox ""
        initGlobal $this-colormap6-listbox ""
        initGlobal $this-colormap1-entry ""
        initGlobal $this-colormap2-entry ""
        initGlobal $this-colormap3-entry ""
        initGlobal $this-colormap4-entry ""
        initGlobal $this-colormap5-entry ""
        initGlobal $this-colormap6-entry ""
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
        global $this-colormap-selection
        global $this-colormap1-listbox
        global $this-colormap2-listbox
        global $this-colormap3-listbox
        global $this-colormap4-listbox
        global $this-colormap5-listbox
        global $this-colormap6-listbox
        global $this-colormap1-entry
        global $this-colormap2-entry
        global $this-colormap3-entry
        global $this-colormap4-entry
        global $this-colormap5-entry
        global $this-colormap6-entry

        sci_toplevel $w 

        wm minsize $w 100 150

        sci_labeledframe $w.frame -labeltext "BUNDLE FIELD OUTPUTS"
        set childframe [$w.frame childsite]
        pack $w.frame -fill both -expand yes

        sci_tabnotebook $childframe.pw  -width 400 -tabpos n
        $childframe.pw add -label "ColorMap1"
        $childframe.pw add -label "ColorMap2" 
        $childframe.pw add -label "ColorMap3" 
        $childframe.pw add -label "ColorMap4"
        $childframe.pw add -label "ColorMap5" 
        $childframe.pw add -label "ColorMap6" 
        $childframe.pw select 0

        pack $childframe.pw -fill both -expand yes

        set colormap1 [$childframe.pw childsite 0]
        set colormap2 [$childframe.pw childsite 1]
        set colormap3 [$childframe.pw childsite 2]
        set colormap4 [$childframe.pw childsite 3]
        set colormap5 [$childframe.pw childsite 4]
        set colormap6 [$childframe.pw childsite 5]

        sci_frame $colormap1.name
        sci_frame $colormap1.sel
        pack $colormap1.name -side top -fill x -expand no -padx 5p
        pack $colormap1.sel -side top -fill both -expand yes -padx 5p

        sci_label $colormap1.name.label -text "Name"
        sci_entry $colormap1.name.entry -textvariable $this-colormap1-name
        set $this-colormap1-entry $colormap1.name.entry
        pack $colormap1.name.label -side left 
        pack $colormap1.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $colormap1.sel.listbox  -selectioncommand [format "%s ChooseColorMap1" $this]
        set $this-colormap1-listbox $colormap1.sel.listbox
        $colormap1.sel.listbox component listbox configure -listvariable $this-colormap-selection -selectmode browse
        pack $colormap1.sel.listbox -fill both -expand yes

        sci_frame $colormap2.name
        sci_frame $colormap2.sel
        pack $colormap2.name -side top -fill x -expand no -padx 5p
        pack $colormap2.sel -side top -fill both -expand yes -padx 5p

        sci_label $colormap2.name.label -text "Name"
        sci_entry $colormap2.name.entry -textvariable $this-colormap2-name
        set $this-colormap2-entry $colormap2.name.entry
        pack $colormap2.name.label -side left 
        pack $colormap2.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $colormap2.sel.listbox  -selectioncommand [format "%s ChooseColorMap2" $this]
        set $this-colormap2-listbox $colormap2.sel.listbox
        $colormap2.sel.listbox component listbox configure -listvariable $this-colormap-selection -selectmode browse
        pack $colormap2.sel.listbox -fill both -expand yes
        
        sci_frame $colormap3.name
        sci_frame $colormap3.sel
        pack $colormap3.name -side top -fill x -expand no -padx 5p
        pack $colormap3.sel -side top -fill both -expand yes -padx 5p

        sci_label $colormap3.name.label -text "Name"
        sci_entry $colormap3.name.entry -textvariable $this-colormap3-name
        set $this-colormap3-entry $colormap3.name.entry
        pack $colormap3.name.label -side left 
        pack $colormap3.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $colormap3.sel.listbox  -selectioncommand [format "%s ChooseColorMap3" $this]
        set $this-colormap3-listbox $colormap3.sel.listbox
        $colormap3.sel.listbox component listbox configure -listvariable $this-colormap-selection -selectmode browse
        pack $colormap3.sel.listbox -fill both -expand yes

        sci_frame $colormap4.name
        sci_frame $colormap4.sel
        pack $colormap4.name -side top -fill x -expand no -padx 5p
        pack $colormap4.sel -side top -fill both -expand yes -padx 5p

        sci_label $colormap4.name.label -text "Name"
        sci_entry $colormap4.name.entry -textvariable $this-colormap4-name
        set $this-colormap4-entry $colormap4.name.entry
        pack $colormap4.name.label -side left 
        pack $colormap4.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $colormap4.sel.listbox  -selectioncommand [format "%s ChooseColorMap4" $this]
        set $this-colormap4-listbox $colormap4.sel.listbox
        $colormap4.sel.listbox component listbox configure -listvariable $this-colormap-selection -selectmode browse
        pack $colormap4.sel.listbox -fill both -expand yes

        sci_frame $colormap5.name
        sci_frame $colormap5.sel
        pack $colormap5.name -side top -fill x -expand no -padx 5p
        pack $colormap5.sel -side top -fill both -expand yes -padx 5p

        sci_label $colormap5.name.label -text "Name"
        sci_entry $colormap5.name.entry -textvariable $this-colormap5-name
        set $this-colormap5-entry $colormap5.name.entry
        pack $colormap5.name.label -side left 
        pack $colormap5.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $colormap5.sel.listbox  -selectioncommand [format "%s ChooseColorMap5" $this]
        set $this-colormap5-listbox $colormap5.sel.listbox
        $colormap5.sel.listbox component listbox configure -listvariable $this-colormap-selection -selectmode browse
        pack $colormap5.sel.listbox -fill both -expand yes
        
        sci_frame $colormap6.name
        sci_frame $colormap6.sel
        pack $colormap6.name -side top -fill x -expand no -padx 5p
        pack $colormap6.sel -side top -fill both -expand yes -padx 5p

        sci_label $colormap6.name.label -text "Name"
        sci_entry $colormap6.name.entry -textvariable $this-colormap6-name
        set $this-colormap6-entry $colormap6.name.entry
        pack $colormap6.name.label -side left 
        pack $colormap6.name.entry -side left -fill x -expand yes
        
        sci_scrolledlistbox $colormap6.sel.listbox  -selectioncommand [format "%s ChooseColorMap6" $this]
        set $this-colormap6-listbox $colormap6.sel.listbox
        $colormap6.sel.listbox component listbox configure -listvariable $this-colormap-selection -selectmode browse
        pack $colormap6.sel.listbox -fill both -expand yes

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
    
    
    method ChooseColorMap1 { } {
        global $this-colormap1-listbox
        global $this-colormap1-name
        global $this-colormap-selection
        
        set colormapnum [[set $this-colormap1-listbox] curselection]
        if [expr [string equal $colormapnum ""] == 0] {
            set $this-colormap1-name  [lindex [set $this-colormap-selection] $colormapnum] 
        }
    }

    method ChooseColorMap2 { } {
        global $this-colormap2-listbox
        global $this-colormap2-name
        global $this-colormap-selection
        
        set colormapnum [[set $this-colormap2-listbox] curselection]
        if [expr [string equal $colormapnum ""] == 0] {
            set $this-colormap2-name  [lindex [set $this-colormap-selection] $colormapnum] 
        }
    }

    method ChooseColorMap3 { } {
        global $this-colormap3-listbox
        global $this-colormap3-name
        global $this-colormap-selection
        
        set colormapnum [[set $this-colormap3-listbox] curselection]
        if [expr [string equal $colormapnum ""] == 0] {
            set $this-colormap3-name  [lindex [set $this-colormap-selection] $colormapnum] 
        }
    }

    method ChooseColorMap4 { } {
        global $this-colormap4-listbox
        global $this-colormap4-name
        global $this-colormap-selection
        
        set colormapnum [[set $this-colormap4-listbox] curselection]
        if [expr [string equal $colormapnum ""] == 0] {
            set $this-colormap4-name  [lindex [set $this-colormap-selection] $colormapnum] 
        }
    }

    method ChooseColorMap5 { } {
        global $this-colormap5-listbox
        global $this-colormap5-name
        global $this-colormap-selection
        
        set colormapnum [[set $this-colormap5-listbox] curselection]
        if [expr [string equal $colormapnum ""] == 0] {
            set $this-colormap5-name  [lindex [set $this-colormap-selection] $colormapnum] 
        }
    }

    method ChooseColorMap6 { } {
        global $this-colormap6-listbox
        global $this-colormap6-name
        global $this-colormap-selection
        
        set colormapnum [[set $this-colormap6-listbox] curselection]
        if [expr [string equal $colormapnum ""] == 0] {
            set $this-colormap6-name  [lindex [set $this-colormap-selection] $colormapnum] 
        }
    }
    
}
