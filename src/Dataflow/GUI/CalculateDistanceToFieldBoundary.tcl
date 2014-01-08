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


itcl::class SCIRun_ChangeFieldData_CalculateDistanceToFieldBoundary {
    inherit Module
    
     constructor { {args ""} } {
        eval configure $args
        set name CalculateDistanceToFieldBoundary
        set_defaults
    }

    method ui {} {
        set w .ui[modname]

        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w

        sci_frame $w.f
        pack $w.f -side top
        sci_checkbutton $w.f.b1 -text "Truncate distance larger than:" -variable $this-truncate
        
        sci_entry $w.f.b3 -textvariable $this-truncate-distance
        
        sci_label $w.f.b4 -text "Data type output field:"
        sci_label $w.f.b6 -text "Basis type output field:"

        grid $w.f.b1 -row 1 -column 0 -sticky w 
        grid $w.f.b3 -row 1 -column 1 -sticky w 
        grid $w.f.b4 -row 2 -column 0 -sticky w
        grid $w.f.b6 -row 3 -column 0 -sticky w

        myselectionbutton $w.f.b5 2 1 { "double" "float" "char" "unsigned char" "short" "unsigned short" "unsigned int" "int"} $this-datatype
        myselectionbutton $w.f.b7 3 1 { "same as input" "constant" "linear" } $this-basistype
        
        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
    
   method myselectionbutton { win x y arglist var} {
        sci_frame $win 
        grid $win  -row $x -column $y -sticky news
        sci_optionmenu $win.c -foreground darkred -command " $this comboget $win.c $var "

        set i 0
        set found 0
        set length [llength $arglist]
        for {set elem [lindex $arglist $i]} {$i<$length} {incr i 1; set elem [lindex $arglist $i]} {
          if {"$elem"=="[set $var]"} {
            set found 1
          }
          $win.c insert end $elem
        }

        if {!$found} {
          $win.c insert end [set $var]
        }

        $win.c select [set $var]  
        pack $win.c	-fill x
    }


    method comboget { win var } {
        if {![winfo exists $win]} {
          return
        }
        if { "$var"!="[$win get]" } {
          set $var [$win get]
        }
    }    
    
}


