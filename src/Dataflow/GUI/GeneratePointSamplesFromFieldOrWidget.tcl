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


itcl::class SCIRun_NewField_GeneratePointSamplesFromFieldOrWidget {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GeneratePointSamplesFromFieldOrWidget
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w

      sci_tabnotebook $w.tabs -raiseselect true \
        -width 350 -height 220 -tabpos n -backdrop gray
      pack $w.tabs -side top -expand 1 -fill both

      set wtab [$w.tabs add -label "Widget" \
          -command "set $this-whichtab Widget"]
      set rtab [$w.tabs add -label "Random" \
          -command "set $this-whichtab Random"]
      if {"[set $this-whichtab]"=="Widget"} {
          $w.tabs view 0
      } else {
          $w.tabs view 1
      }

      sci_labeledframe $wtab.type -labelpos nw \
                       -labeltext "Widget type"
      set type [$wtab.type childsite]
      sci_radiobutton $type.rake -var $this-wtype -value rake -text "Rake" \
          -command "$this-c needexecute" 
      sci_radiobutton $type.ring -var $this-wtype -value ring -text "Ring" \
          -command "$this-c needexecute"
      sci_radiobutton $type.frame -var $this-wtype -value frame -text "Frame" \
          -command "$this-c needexecute"
      pack $type.rake $type.ring $type.frame -side left -padx 5 -pady 5 \
          -fill both -expand yes

      sci_spinner $wtab.f1 -labeltext "Number of samples: " \
          -width 5 -fixed 5 \
	  -validate "$this set-quantity %P $this-widget_seeds" \
          -decrement "$this spin-quantity -1 $wtab.f1 $this-widget_seeds" \
          -increment "$this spin-quantity  1 $wtab.f1 $this-widget_seeds" 

      $wtab.f1 insert 0 [set $this-widget_seeds]

      sci_checkbutton $wtab.auto -text "Execute automatically" \
        -variable $this-autoexecute

      sci_button $wtab.reset -text "Reset Widget" \
          -command "set $this-force-rake-reset 1; $this-c needexecute"

      pack $wtab.type -side top -fill x -pady 5 -anchor w
      pack $wtab.f1 $wtab.auto -side top -pady 5

      pack $wtab.reset -side top -pady 5

      sci_spinner $rtab.f2 -labeltext "Number of samples: " \
        -width 5 -fixed 5 \
        -validate "$this set-quantity %P $this-random_seeds" \
        -decrement "$this spin-quantity -1 $rtab.f2 $this-random_seeds" \
        -increment "$this spin-quantity  1 $rtab.f2 $this-random_seeds" 

      $rtab.f2 insert 0 [set $this-random_seeds]

      pack $rtab.f2 -side top -pady 5

      sci_labeledframe $rtab.dist -labelpos nw \
                       -labeltext "Distribution"
      pack $rtab.dist -fill x -e y
      set dist [$rtab.dist childsite]
      sci_frame $dist.imp 
      sci_frame $dist.uni 
      pack $dist.uni $dist.imp -side left -f both -e y

      sci_label $dist.imp.label -text "Importance Weighted"
      sci_radiobutton $dist.imp.uni -var $this-dist -value impuni \
            -text "Uniform" 
      sci_radiobutton $dist.imp.scat -var $this-dist -value impscat \
            -text "Scattered"
      sci_label $dist.uni.label -text "Not Weighted"
      sci_radiobutton $dist.uni.uni -var $this-dist -value uniuni \
                  -text "Uniform" 
      sci_radiobutton $dist.uni.scat -var $this-dist -value uniscat \
                  -text "Scattered" 
      pack $dist.imp.label $dist.imp.uni $dist.imp.scat \
           $dist.uni.label $dist.uni.uni $dist.uni.scat \
           -side top -padx 5 -pady 2 -anchor w

      sci_checkbutton $rtab.rnginc -text "Increment RNG seed on execute" \
          -var $this-rnginc
      pack $rtab.rnginc -side top -anchor w -padx 8

      sci_frame $rtab.f1 
      pack $rtab.f1 -side top  -anchor w
      sci_label $rtab.f1.rngseed_l -text "Seed value for RNG" -width 23 -anchor w
      sci_entry $rtab.f1.rngseed -text $this-rngseed -width 10
      pack $rtab.f1.rngseed_l $rtab.f1.rngseed -side left -anchor w -padx 8

      sci_frame $rtab.f3
      pack $rtab.f3 -side top -anchor w
      sci_checkbutton $rtab.f3.clamp -text "Clamp to nodes" -var $this-clamp
      pack $rtab.f3.clamp -anchor w -padx 8

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }


    method set-quantity {new quantity} {
      if {! [regexp "\\A\\d*\\.*\\d+\\Z" $new]} {
          return 0
      } elseif {$new < 1} {
          return 0
      } 
      set $quantity $new
      $this-c needexecute
      return 1
    }

    method spin-quantity {step spinner quantity} {
      set newquantity [expr [set $quantity] + $step]

      if {$newquantity < 1} {
          set newquantity 1
      }   
      set $quantity $newquantity
      $spinner delete 0 end
      $spinner insert 0 [set $quantity]
      $this-c needexecute
   }
}
