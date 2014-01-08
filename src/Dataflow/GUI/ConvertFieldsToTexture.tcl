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


catch {rename ConvertFieldsToTexture ""}

itcl::class SCIRun_Visualization_ConvertFieldsToTexture {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
	set name ConvertFieldsToTexture
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w

        set n "$this-c needexecute"
        set s "$this state"

        sci_frame $w.memframe -relief groove -border 2
        sci_label $w.memframe.l -text "Graphics Card Memory (MB)"
        pack $w.memframe -side top -padx 2 -pady 2 -fill both
        pack $w.memframe.l -side top -fill x
        sci_checkbutton $w.memframe.auto -text "Autodetect" -relief flat \
                  -variable $this-card_mem_auto -onvalue 1 -offvalue 0 \
                  -anchor w -command "$s; $n"
        pack $w.memframe.auto -side top -fill x

        sci_frame $w.memframe.bf -relief flat -border 2
              set bf $w.memframe.bf
        pack $bf -side top -fill x
        sci_radiobutton $bf.b0 -text 4 -variable $this-card_mem -value 4 \
            -command $n -state disabled -foreground darkgrey
        sci_radiobutton $bf.b1 -text 8 -variable $this-card_mem -value 8 \
            -command $n -state disabled -foreground darkgrey
        sci_radiobutton $bf.b2 -text 16 -variable $this-card_mem -value 16 \
            -command $n -state disabled -foreground darkgrey
        sci_radiobutton $bf.b3 -text 32 -variable $this-card_mem -value 32 \
            -command $n -state disabled -foreground darkgrey
        sci_radiobutton $bf.b4 -text 64 -variable $this-card_mem -value 64 \
            -command $n -state disabled -foreground darkgrey
        sci_radiobutton $bf.b5 -text 128 -variable $this-card_mem -value 128 \
            -command $n -state disabled -foreground darkgrey
        sci_radiobutton $bf.b6 -text 256 -variable $this-card_mem -value 256 \
            -command $n -state disabled -foreground darkgrey
        sci_radiobutton $bf.b7 -text 512 -variable $this-card_mem -value 512 \
            -command $n -state disabled -foreground darkgrey
        pack $bf.b0 $bf.b1 $bf.b2 $bf.b3 $bf.b4 $bf.b5 $bf.b6 $bf.b7 \
            -side left -expand yes -fill x

        $this state

        sci_frame $w.data
        pack $w.data -side top -fill x

        sci_frame $w.data.scaleframe -relief groove -border 2
        sci_label $w.data.scaleframe.l -text "Scaling"
        pack $w.data.scaleframe -side top -padx 2 -pady 2 -fill both
        pack $w.data.scaleframe.l -side top -fill x

        global $this-is_fixedmin
        sci_frame $w.data.scaleframe.f1 -relief flat
        pack $w.data.scaleframe.f1 -side top -expand yes -fill x
        sci_radiobutton $w.data.scaleframe.f1.b -text "Auto Scale" \
            -variable $this-is_fixed \
            -value 0 -command "$this autoScale"
        pack $w.data.scaleframe.f1.b -side left

        sci_frame $w.data.scaleframe.f2 -relief flat
        pack $w.data.scaleframe.f2 -side top -expand yes -fill x
        sci_radiobutton $w.data.scaleframe.f2.b -text "Fixed Scale" \
          -variable $this-is_fixed \
          -value 1 -command "$this fixedScale"
        pack $w.data.scaleframe.f2.b -side left

        sci_frame $w.data.scaleframe.f3 -relief flat
        pack $w.data.scaleframe.f3 -side top -expand yes -fill x
        sci_label $w.data.scaleframe.f3.l1 -text "value min:  "
        sci_entry $w.data.scaleframe.f3.e1 -textvariable $this-vmin
        sci_label $w.data.scaleframe.f3.l2 -text "value max:  "
        sci_entry $w.data.scaleframe.f3.e2 -textvariable $this-vmax
        pack $w.data.scaleframe.f3.l1 $w.data.scaleframe.f3.e1 \
          $w.data.scaleframe.f3.l2 $w.data.scaleframe.f3.e2 \
          -side left -expand yes -fill x -padx 2 -pady 2

        sci_frame $w.data.scaleframe.f4 -relief flat
        pack $w.data.scaleframe.f4 -side top -expand yes -fill x
        sci_label $w.data.scaleframe.f4.l1 -text " grad min:  "
        sci_entry $w.data.scaleframe.f4.e1 -textvariable $this-gmin
        sci_label $w.data.scaleframe.f4.l2 -text " grad max:  "
        sci_entry $w.data.scaleframe.f4.e2 -textvariable $this-gmax
        pack $w.data.scaleframe.f4.l1 $w.data.scaleframe.f4.e1 \
        $w.data.scaleframe.f4.l2 $w.data.scaleframe.f4.e2 \
        -side left -expand yes -fill x -padx 2 -pady 2

        bind $w.data.scaleframe.f3.e1 <Return> $n
        bind $w.data.scaleframe.f3.e2 <Return> $n

        if { [set $this-is_fixed] } {
            $w.data.scaleframe.f2.b select
            $this fixedScale
        } else {
            $w.data.scaleframe.f1.b select
            $this autoScale
        }

        sci_frame $w.data.histframe -relief groove -border 2
        sci_label $w.data.histframe.l -text "2D Histogram (Multi dimensional transfer function only)"
        pack $w.data.histframe -side top -padx 2 -pady 2 -fill both
        pack $w.data.histframe.l -side top -fill x
        
        sci_frame $w.data.histframe.f1 -relief flat
        pack $w.data.histframe.f1 -side top -expand yes -fill x        
        sci_checkbutton $w.data.histframe.f1.b -text "Create Histogram" \
            -variable $this-histogram
        pack $w.data.histframe.f1.b -side left

        sci_frame $w.data.histframe.f2 -relief flat
        pack $w.data.histframe.f2 -side top -expand yes -fill x        

        sci_label $w.data.histframe.f2.l1 -text "histogram gamma"
        sci_entry $w.data.histframe.f2.e1 -textvariable $this-gamma
        pack $w.data.histframe.f2.l1 $w.data.histframe.f2.e1  -side left \


        makeSciButtonPanel $w $w $this
        moveToCursor $w

    }

    method autoScale { } {
        global $this-is_fixed
        set w .ui[modname]
        
        set $this-is_fixed 0

        set color "#505050"

        $w.data.scaleframe.f3.l1 configure -foreground $color
        $w.data.scaleframe.f3.e1 configure -state disabled -foreground $color
        $w.data.scaleframe.f3.l2 configure -foreground $color
        $w.data.scaleframe.f3.e2 configure -state disabled -foreground $color

        $w.data.scaleframe.f4.l1 configure -foreground $color
        $w.data.scaleframe.f4.e1 configure -state disabled -foreground $color
        $w.data.scaleframe.f4.l2 configure -foreground $color
        $w.data.scaleframe.f4.e2 configure -state disabled -foreground $color
   }	

    method fixedScale { } {
        global $this-is_fixed
        set w .ui[modname]

        set $this-is_fixed 1


        $w.data.scaleframe.f3.l1 configure -foreground black
        $w.data.scaleframe.f3.e1 configure -state normal -foreground black
        $w.data.scaleframe.f3.l2 configure -foreground black
        $w.data.scaleframe.f3.e2 configure -state normal -foreground black

        $w.data.scaleframe.f4.l1 configure -foreground black
        $w.data.scaleframe.f4.e1 configure -state normal -foreground black
        $w.data.scaleframe.f4.l2 configure -foreground black
        $w.data.scaleframe.f4.e2 configure -state normal -foreground black        
    }

    method set_card_mem {mem} {
        set $this-card_mem $mem
    }

    method set_card_mem_auto {auto} {
        set $this-card_mem_auto $auto
    }

    method state {} {
      set w .ui[modname]
      if {[winfo exists $w] == 0} {
          return
      }
      if {[set $this-card_mem_auto] == 1} {
                $this deactivate $w.memframe.bf.b0
                $this deactivate $w.memframe.bf.b1
                $this deactivate $w.memframe.bf.b2
                $this deactivate $w.memframe.bf.b3
                $this deactivate $w.memframe.bf.b4
                $this deactivate $w.memframe.bf.b5
                $this deactivate $w.memframe.bf.b6
                $this deactivate $w.memframe.bf.b7
      } else {
                $this activate $w.memframe.bf.b0
                $this activate $w.memframe.bf.b1
                $this activate $w.memframe.bf.b2
                $this activate $w.memframe.bf.b3
                $this activate $w.memframe.bf.b4
                $this activate $w.memframe.bf.b5
                $this activate $w.memframe.bf.b6
                $this activate $w.memframe.bf.b7
      }
    }

    method activate { w } {
      $w configure -state normal -foreground black
    }

    method deactivate { w } {
      $w configure -state disabled -foreground darkgrey
    }

}
