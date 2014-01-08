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
#    File   : SetRegionAttribs.tcl
#    Author : Martin Cole
#    Date   : Tue Mar 28 14:00:05 2006

catch {rename SCIRun_NewField_InterfaceWithTetGen ""}

itcl::class SCIRun_NewField_InterfaceWithTetGen {
  inherit Module

   constructor { {args ""} } {
        eval configure $args
    set name InterfaceWithTetGen
    global $this-minRadius
    global $this-maxVolConstraint
  }

  method ui {} {
    set w .ui[modname]
    if {[winfo exists $w]} {
      raise $w
      return;
    }


    sci_toplevel $w
    wm minsize $w 400 280
    sci_label $w.l -text "TetGen Options"
    sci_frame $w.switches -relief groove -borderwidth 2

    sci_frame $w.switches.pOption -relief flat
    # option -p:
    sci_checkbutton $w.switches.pOption.piecewise \
      -text "-p: Tetrahedralize a piecewise linear complex (PLC)" \
      -onvalue 1 -offvalue 0 -variable $this-piecewiseFlag
    pack $w.switches.pOption.piecewise -side top -anchor nw

    sci_frame $w.switches.attribOption -relief flat
    # option -A:
    sci_checkbutton $w.switches.attribOption.attribs \
      -text "-A: Assign attributes to identify tetrahedra in certain regions" \
      -onvalue 1 -offvalue 0 -variable $this-assignFlag \
      -command "$this toggleAttribOption"
    sci_frame $w.switches.attribOption.subOption -relief flat -borderwidth 2
    sci_checkbutton $w.switches.attribOption.subOption.nonzeroAttrib \
      -text "-AA: Assign non-zero attributes" \
      -onvalue 1 -offvalue 0 -variable $this-setNonzeroAttributeFlag \
      -state disabled
    pack $w.switches.attribOption.attribs -side top -anchor nw
    pack $w.switches.attribOption.subOption.nonzeroAttrib -anchor nw
    pack $w.switches.attribOption.subOption -anchor s
     
    sci_frame $w.switches.yOption -relief flat
    # option -Y:
    sci_checkbutton $w.switches.yOption.suppressBoundSplit \
      -text "-Y: Suppress boundary facets/segments splitting" \
      -onvalue 1 -offvalue 0 -variable $this-suppressSplitFlag \
      -command "$this toggleSuppressOption"
    sci_frame $w.switches.yOption.subOption -relief flat -borderwidth 2
    sci_checkbutton $w.switches.yOption.subOption.suppressAllSplit \
      -text "-YY: Suppress all boundary splitting" \
      -onvalue 1 -offvalue 0 -variable $this-setSplitFlag \
      -state disabled
    pack $w.switches.yOption.subOption.suppressAllSplit -side top -anchor nw
    pack $w.switches.yOption.suppressBoundSplit -anchor nw
    pack $w.switches.yOption.subOption -anchor s
    
    sci_frame $w.switches.zOption -relief flat
    # option -z:
    sci_checkbutton $w.switches.zOption.numberOutput \
      -text "-z: Number all output items starting from zero" \
      -onvalue 1 -offvalue 0 -variable $this-numberOutputFlag
    pack $w.switches.zOption.numberOutput -side top -anchor nw
     
    sci_frame $w.switches.qOption -relief flat
    # option -q:
    sci_checkbutton $w.switches.qOption.quality \
      -text "-q: Quality mesh generation" \
      -onvalue 1 -offvalue 0 -variable $this-qualityFlag \
      -command "$this toggleQualityOption"
    sci_frame $w.switches.qOption.subOption -relief flat -borderwidth 2
    sci_checkbutton $w.switches.qOption.subOption.setRatio \
      -text "Specify minimum radius-edge ratio" \
      -onvalue 1 -offvalue 0 -variable $this-setRatioFlag \
      -state disabled \
      -command "$this toggleRatioBound"
    sci_entry $w.switches.qOption.subOption.minRadiusEdgee -width 10 \
      -state disabled -textvariable $this-minRadius
    pack $w.switches.qOption.subOption.setRatio -side top -anchor nw
    pack $w.switches.qOption.subOption.minRadiusEdgee -side bottom
    pack $w.switches.qOption.quality -anchor nw
    pack $w.switches.qOption.subOption -anchor s
     
    sci_frame $w.switches.aOption -relief flat
    # option -a:
    sci_checkbutton $w.switches.aOption.constraint \
      -text "-a: Impose volume constraint" \
      -onvalue 1 -offvalue 0 -variable $this-volConstraintFlag \
      -command "$this toggleVolumeConstraint"
    sci_frame $w.switches.aOption.subOption -relief flat -borderwidth 2
    sci_checkbutton $w.switches.aOption.subOption.maxVolConstraint \
      -text "Specify a maximum volume constraint on all tetrahedra" \
      -onvalue 1 -offvalue 0 -variable $this-setMaxVolConstraintFlag \
      -state disabled \
      -command "$this toggleAssignAttributes"
    sci_entry $w.switches.aOption.subOption.maxVolConstrainte -width 10 \
      -state disabled -textvariable $this-maxVolConstraint
    pack $w.switches.aOption.subOption.maxVolConstrainte -side bottom
    pack $w.switches.aOption.constraint -anchor nw
    pack $w.switches.aOption.subOption.maxVolConstraint
    pack $w.switches.aOption.subOption -anchor s
    
    sci_frame $w.switches.dOption -relief flat
    # option -d:
    sci_checkbutton $w.switches.dOption.detectIntersections \
      -text "-d: Detect intersections of PLC facets" \
      -onvalue 1 -offvalue 0 -variable $this-detectIntersectionsFlag
    pack $w.switches.dOption.detectIntersections -anchor nw
    
    sci_frame $w.moreSwitches -relief groove -borderwidth 2
    sci_label $w.moreSwitches.l \
      -text "Additional Tetgen command line options"
    sci_entry $w.moreSwitches.moreSwitchese -width 20 \
      -state normal -textvariable $this-moreSwitches
    pack $w.moreSwitches.l -side top -fill x -expand yes
    pack $w.moreSwitches.moreSwitchese -anchor s
     
    pack $w.l -side top
    pack $w.switches.pOption \
         $w.switches.attribOption \
         $w.switches.zOption \
         $w.switches.yOption \
         $w.switches.qOption \
         $w.switches.aOption \
         $w.switches.dOption \
         -fill x -anchor nw -expand no
    
    pack $w.switches -side top -fill x -expand yes
    pack $w.moreSwitches -side top -fill x -expand yes
    
    makeSciButtonPanel $w $w $this
    toggleQualityOption
    toggleAttribOption
    toggleSuppressOption
    toggleVolumeConstraint
    moveToCursor $w
  }

  method toggleQualityOption {} {
    set w .ui[modname]
    if { [set $this-qualityFlag] == 1 } {
      $w.switches.qOption.subOption.setRatio configure -state normal
    } else {
      set $this-setRatioFlag 0
      $w.switches.qOption.subOption.setRatio configure -state disabled
      toggleRatioBound
    }
  }

  method toggleAttribOption {} {
    set w .ui[modname]
    if { [set $this-assignFlag] == 1 } {
      $w.switches.attribOption.subOption.nonzeroAttrib configure -state normal
    } else {
      set $this-setNonzeroAttributeFlag 0
      $w.switches.attribOption.subOption.nonzeroAttrib configure -state disabled
    }
  }


  method toggleSuppressOption {} {
    set w .ui[modname]
    if { [set $this-suppressSplitFlag] == 1 } {
      $w.switches.yOption.subOption.suppressAllSplit configure -state normal
    } else {
      set $this-setSplitFlag 0
      $w.switches.yOption.subOption.suppressAllSplit configure -state disabled
    }
  }
  
  method toggleVolumeConstraint {} {
    set w .ui[modname]
    if { [set $this-volConstraintFlag] == 1 } {
      $w.switches.aOption.subOption.maxVolConstraint configure -state normal
    } else {
      set $this-setMaxVolConstraintFlag 0
      $w.switches.aOption.subOption.maxVolConstraint configure -state disabled
      toggleAssignAttributes
    }
  }

  method toggleRatioBound {} {
    set w .ui[modname]
    if { [set $this-setRatioFlag] == 1 } {
      $w.switches.qOption.subOption.minRadiusEdgee configure -state normal
    } else {
      $w.switches.qOption.subOption.minRadiusEdgee configure -state disabled
    }
  }

  method toggleAssignAttributes {} {
    set w .ui[modname]
    if { [set $this-setMaxVolConstraintFlag] == 1 } {
      $w.switches.aOption.subOption.maxVolConstrainte configure -state normal
    } else {
      $w.switches.aOption.subOption.maxVolConstrainte configure -state disabled
    }
  }
}
