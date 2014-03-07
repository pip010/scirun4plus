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

itcl::class SCIRun_ChangeFieldData_CalculateIsInsideField {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name CalculateIsInsideField
        set_defaults
    }
    
    method set_defaults {} {
      global $this-outputtype

      global $this-outval
      global $this-inval
      
      set $this-outputtype "double"
      set $this-outval 0
      set $this-inval  1
      
    }
    
    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        
        sci_toplevel $w
        
        sci_frame $w.f
        pack $w.f
        sci_label $w.f.lab1 -text "Data location"
        grid $w.f.lab1 -row 0 -column 0 -sticky e
        sci_label $w.f.lab2 -text "Data type"
        grid $w.f.lab2 -row 1 -column 0 -sticky e
        sci_label $w.f.lab3 -text "Outside value"
        grid $w.f.lab3 -row 2 -column 0 -sticky e
        sci_label $w.f.lab4 -text "Inside value"
        grid $w.f.lab4 -row 3 -column 0 -sticky e
        
        
        myselectionbutton $w.f.sel2 1 1 { "same as input" "char" "short" "unsigned short" "unsigned int" "int" "float" "double" } $this-outputtype

        sci_entry $w.f.e1 -textvariable $this-outval
        sci_entry $w.f.e2 -textvariable $this-inval        
        grid $w.f.e1 -row 2 -column 1 -sticky news
        grid $w.f.e2 -row 3 -column 1 -sticky news
        
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

    
