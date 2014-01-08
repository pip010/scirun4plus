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


#  ColorPicker.tcl
#  Written by:
#   James Purciful
#   Department of Computer Science
#   University of Utah
#   Mar. 1994


# If the window passed in ('w') does not exist, makeColorPicker will
# create the window for you.

proc makeColorPicker {w var command cancel} {
  global $var
  global Color
  global $var-r $var-g $var-b $var-a
  global $w-r $w-g $w-b $w-a
  global $w-rgbhsv

  if ![winfo exists $w] {
    sci_toplevel $w
    wm title $w "Color Picker" 
    # Withdraw the window until it has been completely created.
    # This will keep it from flickering when it is repositioned
    # (to the cursor).
    wm withdraw $w 
  }

  set $w-rgbhsv "rgb"

  set $w-r [set $var-r]
  set $w-g [set $var-g]
  set $w-b [set $var-b]
  if [info exists $var-a] {
    set $w-a [set $var-a]
  }

  sci_frame $w.c

  set ir [expr int([set $w-r] * 65535)]
  set ig [expr int([set $w-g] * 65535)]
  set ib [expr int([set $w-b] * 65535)]
  
  
  sci_frame $w.c.picks
  set picks $w.c.picks  
  sci_frame $w.c.f 
  sci_frame $w.c.f.opts2 
  sci_frame $w.c.f.col -relief ridge -borderwidth 4 -height 1.5c -width 4c \
    -background [format #%04x%04x%04x $ir $ig $ib]

  set col $w.c.f.col

  sci_frame $picks.rgb -relief groove -borderwidth 0
  set rgb $picks.rgb
  sci_scale $rgb.s1 -label "Red" -from 0.0 -to 1.0 -length 6c -showvalue true \
    -orient horizontal -resolution .01 \
    -digits 3 -variable $w-r 
  sci_scale $rgb.s2 -label "Green" -from 0.0 -to 1.0 -length 6c -showvalue true \
    -orient horizontal -resolution .01 \
    -digits 3 -variable $w-g 
  sci_scale $rgb.s3 -label "Blue" -from 0.0 -to 1.0 -length 6c -showvalue true \
    -orient horizontal -resolution .01 \
    -digits 3 -variable $w-b
  pack $rgb.s1 -side top -padx 2 -pady 2 -anchor nw -fill y
  pack $rgb.s2 -side top -padx 2 -pady 2 -anchor nw -fill y
  pack $rgb.s3 -side top -padx 2 -pady 2 -anchor nw -fill y
  $rgb.s1 set [set $w-r]
  $rgb.s2 set [set $w-g]
  $rgb.s3 set [set $w-b]

  if [info exists $var-a] {
    sci_scale $rgb.s4 -label Alpha -from 0.0 -to 1.0 -length 6c \
      -showvalue true -orient horizontal -resolution .01 \
      -digits 3 -variable $w-a
    pack $rgb.s4 -in $picks.rgb -side top -padx 2 -pady 2 \
      -anchor nw -fill y
    $rgb.s4 set [set $w-a]
  }

  frame $picks.hsv -background $Color(MenuBackGround)
  set hsv $picks.hsv
  sci_scale $hsv.s1 -label Hue -from 0.0 -to 360.0 -length 6c -showvalue true \
    -orient horizontal -resolution .01 \
    -digits 3 -variable $w-h
  sci_scale $hsv.s2 -label Saturation -from 0.0 -to 1.0 -length 6c -showvalue true \
    -orient horizontal -resolution .01 \
    -digits 3 -variable $w-s
  sci_scale $hsv.s3 -label Value -from 0.0 -to 1.0 -length 6c -showvalue true \
    -orient horizontal -resolution .01 \
    -digits 3 -variable $w-v)
  pack $hsv.s1 -side top -padx 2 -pady 2 -anchor nw -fill y
  pack $hsv.s2 -side top -padx 2 -pady 2 -anchor nw -fill y
  pack $hsv.s3 -side top -padx 2 -pady 2 -anchor nw -fill y

  if [info exists $var-a] {
    sci_scale $hsv.s4 -label Alpha -from 0.0 -to 1.0 -length 6c \
      -showvalue true -orient horizontal -resolution .01 \
      -digits 3 -variable $w-a
    pack $hsv.s4 -in $picks.hsv -side top -padx 2 -pady 2 -anchor nw \
      -fill y
  }

  $rgb.s1 configure -command "cpsetrgb $col $rgb.s1 $rgb.s2 $rgb.s3 \
  $hsv.s1 $hsv.s2 $hsv.s3 "
  $rgb.s2 configure -command "cpsetrgb $col $rgb.s1 $rgb.s2 $rgb.s3 \
  $hsv.s1 $hsv.s2 $hsv.s3 "
  $rgb.s3 configure -command "cpsetrgb $col $rgb.s1 $rgb.s2 $rgb.s3 \
  $hsv.s1 $hsv.s2 $hsv.s3 "
  $hsv.s1 configure -command "cpsethsv $col $rgb.s1 $rgb.s2 $rgb.s3 \
  $hsv.s1 $hsv.s2 $hsv.s3 "
  $hsv.s2 configure -command "cpsethsv $col $rgb.s1 $rgb.s2 $rgb.s3 \
  $hsv.s1 $hsv.s2 $hsv.s3 "
  $hsv.s3 configure -command "cpsethsv $col $rgb.s1 $rgb.s2 $rgb.s3 \
  $hsv.s1 $hsv.s2 $hsv.s3 "

  sci_radiobutton $w.c.f.opts2.rgb -text RGB -variable $w-rgbhsv -value rgb \
    -command "cptogrgbhsv $w $picks $rgb $hsv" -background $Color(MenuBackGround)
  sci_radiobutton $w.c.f.opts2.hsv -text HSV -variable $w-rgbhsv -value hsv \
    -command "cptogrgbhsv $w $picks $rgb $hsv" -background $Color(MenuBackGround)

  pack $w.c.f.opts2.hsv $w.c.f.opts2.rgb -side top -padx 2 -pady 2 -anchor e
  pack $w.c.f.opts2 -side right -padx 2 -anchor e
  pack $w.c.f.col -side right -fill x -padx 2 

  frame $w.c.opts -background $Color(MainBackGround) -bd 1
  
  sci_button $w.c.opts.apply -text "Apply" \
    -command "cpcommitcolor $var $rgb.s1 $rgb.s2 $rgb.s3 $rgb.s4 \"$command\""
  
  sci_button $w.c.opts.ok -text "OK" \
    -command "cpcommitcolor $var $rgb.s1 $rgb.s2 $rgb.s3 $rgb.s4 \"$command\"; $cancel"
  
  sci_button $w.c.opts.cancel -text "Cancel" -command $cancel
  
  pack $w.c.opts.cancel  $w.c.opts.apply $w.c.opts.ok -side left -anchor e -padx 2 -pady 2

  if { [set $w-rgbhsv] == "rgb" } {
    pack $rgb -in $picks -side left -padx 2 -pady 2 -expand 1 -fill x
  }
  if { [set $w-rgbhsv] == "hsv" } {
    pack $hsv -in $picks -side left -padx 2 -pady 2 -expand 1 -fill x
  }

  pack $picks $w.c.f -side top \
    -padx 2 -pady 5 -expand 1 -fill both
  pack $w.c.opts -side top -expand 1 -fill x -anchor e 
  pack $w.c

  TooltipMultiWidget "$w.c.f.opts2.rgb $w.c.f.opts2.hsv" "These radio buttons allow you to select Red/Blue/Green (RGB)\nor Hue/Saturation/Value (HSV) methods of entering a color.\nTo select a shade of gray, it is easier to use the HSV mode."

  if { [winfo toplevel $w] == $w } {
      # If $w is a toplevel window, move it to the cursor.
    wm deiconify $w
  }
}

proc Max {n1 n2 n3} {
  if [expr $n1 >= $n2] {
    if [expr $n1 >= $n3] {
      return $n1
    } else {
      return $n3
    }
  } else {
    if [expr $n2 >= $n3] {
      return $n2
    } else {
      return $n3
    }
  }
}

proc Min {n1 n2 n3} {
  if [expr $n1 <= $n2] {
    if [expr $n1 <= $n3] {
	    return $n1
    } else {
	    return $n3
    }
  } else {
    if [expr $n2 <= $n3] {
	    return $n2
    } else {
	    return $n3
    }
  }
}

proc cpsetrgb {col rs gs bs hs ss vs val} {
  # Do inverse transformation to HSV
  set max [Max [$rs get] [$gs get] [$bs get]]
  set min [Min [$rs get] [$gs get] [$bs get]]
  # $ss set [expr ($max == 0.0) ? 0.0 : (($max-$min)/$max)]
  if {$max == 0.0} {
    $ss set 0.0
  } else {
    $ss set [expr ($max-$min)/$max]
  }
  if [expr [$ss get] != 0.0] {
    set rl [expr ($max-[$rs get])/($max-$min)]
    set gl [expr ($max-[$gs get])/($max-$min)]
    set bl [expr ($max-[$bs get])/($max-$min)]
    if [expr $max == [$rs get]] {
	    if [expr $min == [$gs get]] {
        $hs set [expr 60.0*(5.0+$bl)]
	    } else {
        $hs set [expr 60.0*(1.0-$gl)]
	    }
    } elseif [expr $max == [$gs get]] {
	    if [expr $min == [$bs get]] {
        $hs set [expr 60.0*(1.0+$rl)]
	    } else {
        $hs set [expr 60.0*(3.0-$bl)]
	    }
    } else {
	    if [expr $min == [$rs get]] {
        $hs set [expr 60.0*(3.0+$gl)]
	    } else {
        $hs set [expr 60.0*(5.0-$rl)]
	    }
    }
  } else {
    $hs set 0.0
  }
  $vs set $max

  cpsetcol $col [$rs get] [$gs get] [$bs get]
  update idletasks
}

proc cpsethsv {col rs gs bs hs ss vs val} {
  # Convert to RGB...
  while {[$hs get] >= 360.0} {
    $hs set [expr [$hs get] - 360.0]
  }
  while {[$hs get] < 0.0} {
    $hs set [expr [$hs get] + 360.0]
  }
  set h6 [expr [$hs get]/60.0]
  set i [expr int($h6)]
  set f [expr $h6-$i]
  set p1 [expr [$vs get]*(1.0-[$ss get])]
  set p2 [expr [$vs get]*(1.0-([$ss get]*$f))]
  set p3 [expr [$vs get]*(1.0-([$ss get]*(1-$f)))]
  switch $i {
    0 {$rs set [$vs get] ; $gs set $p3 ; $bs set $p1}
    1 {$rs set $p2 ; $gs set [$vs get] ; $bs set $p1}
    2 {$rs set $p1 ; $gs set [$vs get] ; $bs set $p3}
    3 {$rs set $p1 ; $gs set $p2 ; $bs set [$vs get]}
    4 {$rs set $p3 ; $gs set $p1 ; $bs set [$vs get]}
    5 {$rs set [$vs get] ; $gs set $p1 ; $bs set $p2}
    default {$rs set 0 ; $gs set 0 ; $bs set 0}
  }

  cpsetcol $col [$rs get] [$gs get] [$bs get]
  update idletasks
}

proc cpsetcol {col r g b} {
  set ir [expr int($r * 65535)]
  set ig [expr int($g * 65535)]
  set ib [expr int($b * 65535)]

  $col config -background [format #%04x%04x%04x $ir $ig $ib]
}

proc cpcommitcolor {var rs gs bs as command} {
  global $var-r $var-g $var-b $var-a
  set $var-r [$rs get]
  set $var-g [$gs get]
  set $var-b [$bs get]
  if [info exists $var-a] {
    set $var-a [$as get]
  }
  eval $command
}

proc cptogrgbhsv {w picks rgb hsv} {
  global $w-rgbhsv
  if { [set $w-rgbhsv] == "rgb" } {
    pack forget $hsv
    pack $rgb -in $picks -side left
  } else {
    pack forget $rgb
    pack $hsv -in $picks -side left
  }
}
