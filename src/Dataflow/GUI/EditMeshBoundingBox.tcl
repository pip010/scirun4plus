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


itcl::class SCIRun_ChangeMesh_EditMeshBoundingBox {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name EditMeshBoundingBox

	# The width of the first column of the data display.
	setGlobal $this-firstwidth 12
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w
        sci_labeledframe $w.att -labelpos nw \
                         -labeltext "Input Field Attributes" 
                   
        pack $w.att -side top -fill x -expand y
        set att [$w.att childsite]
	
        labelpair3 $att.l1 "Center (x,y,z)" \
          $this-inputcenterx $this-inputcentery $this-inputcenterz
        labelpair3 $att.l2 "Size (x,y,z)" \
          $this-inputsizex $this-inputsizey $this-inputsizez
        pack $att.l1 $att.l2 -side top -fill x

        sci_labeledframe $w.edit -labelpos nw \
                         -labeltext "Output Field Attributes" 
        pack $w.edit -side top
        set edit [$w.edit childsite]
	
        labelentry3 $edit.l1 "Center (x,y,z)" \
          $this-outputcenterx $this-outputcentery $this-outputcenterz \
          "$this-c needexecute" \
          $this-useoutputcenter
        labelentry3 $edit.l2 "Size (x,y,z)" \
          $this-outputsizex $this-outputsizey \
          $this-outputsizez "$this-c needexecute" \
          $this-useoutputsize

        pack $edit.l1 $edit.l2 -side top 


        sci_labeledframe $w.scale -labelpos nw -labeltext "Widget Scale/Mode" 
        set scale [$w.scale childsite]
        sci_frame $scale.f 
        pack $scale.f -side left
        
        sci_label  $scale.f.l1 -text "SCALE:"
        grid   $scale.f.l1 -column 0 -row 0
        sci_frame  $scale.f.f1
        grid   $scale.f.f1 -column 1 -row 0
        sci_button $scale.f.f1.incr -text "++" -command "$this-c scale 1.25"
        sci_button $scale.f.f1.incr2 -text "+" -command "$this-c scale 1.05"
        sci_button $scale.f.f1.decr -text "-" -command "$this-c scale [expr 1.0/1.05]"
        sci_button $scale.f.f1.decr2 -text "--" -command "$this-c scale [expr 1.0/1.25]"
        pack   $scale.f.f1.incr $scale.f.f1.incr2 $scale.f.f1.decr $scale.f.f1.decr2 -side left -anchor w

        sci_label  $scale.f.l2 -text "MODE:"
        grid   $scale.f.l2 -column 2 -row 0
        sci_button $scale.f.nextmode -text "NextMode" -command "$this-c nextmode"
        grid   $scale.f.nextmode -column 3 -row 0

        sci_radiobutton $scale.f.r1 -variable $this-restrict-translation -value 0 -text "NO TRANSLATION RESTRICTION" -command "$this-c restricttranslation 0"
        grid $scale.f.r1 -column 0 -row 1 -columnspan 4 -sticky w
        sci_radiobutton $scale.f.r2 -variable $this-restrict-translation -value 1 -text "XYZ TRANSLATION RESTRICTION" -command "$this-c restricttranslation 1"
        grid $scale.f.r2 -column 0 -row 2 -columnspan 4 -sticky w
        sci_radiobutton $scale.f.r3 -variable $this-restrict-translation -value 2 -text "RDI TRANSLATION RESTRICTION" -command "$this-c restricttranslation 2"
        grid $scale.f.r3 -column 0 -row 3 -columnspan 4 -sticky w

        sci_checkbutton $scale.f.rx -variable $this-restrict-x -command "$this-c restrictx" -text "X"
        sci_checkbutton $scale.f.ry -variable $this-restrict-y -command "$this-c restricty" -text "Y"
        sci_checkbutton $scale.f.rz -variable $this-restrict-z -command "$this-c restrictz" -text "Z"
        grid $scale.f.rx -column 4 -row 2
        grid $scale.f.ry -column 5 -row 2
        grid $scale.f.rz -column 6 -row 2

        sci_checkbutton $scale.f.rr -variable $this-restrict-r -command "$this-c restrictr" -text "R"
        sci_checkbutton $scale.f.rd -variable $this-restrict-d -command "$this-c restrictd" -text "D"
        sci_checkbutton $scale.f.ri -variable $this-restrict-i -command "$this-c restricti" -text "I"
        grid $scale.f.rr -column 4 -row 3
        grid $scale.f.rd -column 5 -row 3
        grid $scale.f.ri -column 6 -row 3

        pack $w.scale -side top -fill x -expand y

        makeSciButtonPanel $w $w $this \
            "\"Reset Widget\" \"$this reset\" \"\"" \
            "\"In to Out\" \"$this copy_attributes\" \"Copies the Input Field Attribute values\nto the Output Field Attribute text fields.\n(This is just for user convenience.)\" "
        moveToCursor $w	
    }

    method reset {} {
        global $this-resetting
        set $this-resetting 1
        $this-c needexecute
    }

    method copy_attributes {} {
        set w .ui[modname]
        if {![winfo exists $w]} {
            return
        }
        set att [$w.att childsite]
        set edit [$w.edit childsite]
        set $this-outputcenterx [set $this-inputcenterx]
        set $this-outputcentery [set $this-inputcentery]
        set $this-outputcenterz [set $this-inputcenterz]
        set $this-outputsizex [set $this-inputsizex]
        set $this-outputsizey [set $this-inputsizey]
        set $this-outputsizez [set $this-inputsizez]
    }
}
