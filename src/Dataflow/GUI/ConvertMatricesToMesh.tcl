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


itcl::class SCIRun_NewField_ConvertMatricesToMesh {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name ConvertMatricesToMesh

	# The width of the first column of the data display.
	setGlobal $this-firstwidth 12
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w
        wm maxsize $w 700 215 

        sci_labeledframe $w.att -labelpos nw \
                         -labeltext "Input Field Type" 
                   
        pack $w.att -anchor nw
        set att [$w.att childsite]
        
        labelpair $att.l1 "Name" $this-fldname
        labelpair $att.l2 "Typename" $this-inputdatatype
        pack $att.l1 $att.l2 -side top -anchor nw

        sci_labeledframe $w.fbt -labelpos nw \
                         -labeltext "Output Mesh Type"
        pack $w.fbt -anchor nw
        set fbt [$w.fbt childsite]
        labelcombo $fbt.l1 "Mesh Type" \
            { \
            Curve \
            HexVol \
            PointCloud \
            PrismVol \
            QuadSurf \
            TetVol \
            TriSurf \
              } \
            $this-fieldbasetype

        pack $fbt.l1 -side top -anchor nw

        sci_labeledframe $w.edit -labelpos nw \
                         -labeltext "Output Field Type" 
        pack $w.edit -anchor nw
        set edit [$w.edit childsite]
        labelcombo $edit.l1 "Data Type" \
          {"unsigned char" "unsigned short" "unsigned int" \
          "char" "short" "int" "float" "double" "Vector" "Tensor"} \
             $this-datatype
        pack $edit.l1 -side top -anchor nw

        makeSciButtonPanel $w $w $this
        moveToCursor $w
      }

      method labelcombo { win text1 arglist var} {
        sci_frame $win 
        pack $win -side top -padx 5 -anchor nw
        sci_label $win.l1 -text $text1 -width [set $this-firstwidth] \
                -anchor w -just left
        sci_label $win.colon  -text ":" -width 2 -anchor w -just left
        sci_optionmenu $win.c -foreground darkred \
          -command " $this comboget $win.c $var "

        set i 0
        set found 0
        set length [llength $arglist]
        for {set elem [lindex $arglist $i]} {$i<$length} \
            {incr i 1; set elem [lindex $arglist $i]} {
            if {"$elem"=="[set $var]"} {
          set found 1
            }
            $win.c insert end $elem
        }

        if {!$found} {
            $win.c insert end [set $var]
        }

        sci_label $win.l2 -text "" -width 40 -anchor w -just left

        # hack to associate optionmenus with a textvariable
        bind $win.c <Map> "$win.c select {[set $var]}"

        pack $win.l1 $win.colon -side left -anchor nw
        pack $win.c $win.l2 -side left  -anchor nw
     }

     method comboget { win var } {
        if {![winfo exists $win]} {
            return
        }
        if { "$var"!="[$win get]" } {
            set $var [$win get]
        }
    }

    method config_labelcombo { win arglist sel} {
	if {![winfo exists $win]} {
	    return
	}
	$win.c delete 0 end
	if {[llength $arglist]==0} {
	    $win.c insert end ""
	}
	set i 0
	set length [llength $arglist]
	for {set elem [lindex $arglist $i]} {$i<$length} \
	    {incr i 1; set elem [lindex $arglist $i]} {
	    $win.c insert end $elem
	}
	
	if {"$sel"!="---"} {
	    $win.c select $sel
	}
    }
}




