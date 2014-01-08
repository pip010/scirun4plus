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

#
#  CreateMeshFromNrrd.tcl:
#
#  Create a Field Mesh from Nrrd Daata. Incoming Nrrds may consist of
#  mesh points and optionally mesh connections.
#
#  Written by:
#   Allen R. Sanderson
#   SCI Intitute
#   University of Utah
#   Febuary 2007
#

catch {rename Teem_Converters_CreateMeshFromNrrd ""}

itcl::class Teem_Converters_CreateMeshFromNrrd {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name CreateMeshFromNrrd
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        sci_toplevel $w

        sci_frame $w.f
        pack $w.f -padx 2 -pady 2 -side top -expand yes
        
        sci_frame $w.f.options
        pack $w.f.options -side top -expand yes

        sci_labeledframe $w.f.options.quadtet \
            -labelpos nw -labeltext "Unstructured Cell Type when\nPoints per Connection = 4:"
        pack $w.f.options.quadtet -side top -expand yes -fill x

        set quadtet [$w.f.options.quadtet childsite]

        sci_radiobutton $quadtet.auto -text "Auto" \
            -variable $this-quad-or-tet -value "Auto"
        sci_radiobutton $quadtet.tet -text "Tet" \
            -variable $this-quad-or-tet -value "Tet"
        sci_radiobutton $quadtet.quad -text "Quad" \
            -variable $this-quad-or-tet -value "Quad"

        pack $quadtet.auto $quadtet.tet $quadtet.quad \
            -side left -anchor nw -padx 3


        sci_labeledframe $w.f.options.pccurve \
            -labelpos nw -labeltext "Structured/Unstructured Ambiguity:"
        pack $w.f.options.pccurve -side top -expand yes -fill x
        set pccurve [$w.f.options.pccurve childsite]

        sci_radiobutton $pccurve.auto -text "Auto" \
            -variable $this-struct-or-unstruct -value "Auto"
        sci_radiobutton $pccurve.pc -text "Point Cloud" \
            -variable $this-struct-or-unstruct -value "PointCloud"
        sci_radiobutton $pccurve.curve -text "Struct Curve" \
            -variable $this-struct-or-unstruct -value "StructCurve"
        pack $pccurve.auto $pccurve.pc $pccurve.curve \
            -side left -anchor nw -padx 3

        # Input Dataset
        sci_label $w.f.options.datasetslab -text "Datasets:"
        pack $w.f.options.datasetslab -side top -anchor nw -pady 3

        sci_frame $w.f.options.datasets
        pack $w.f.options.datasets -side top -anchor nw -pady 5

        global $this-datasets
        set_names [set $this-datasets]

        makeSciButtonPanel $w $w $this
        moveToCursor $w

        pack $w.f -expand 1 -fill x
    }

    method set_names {datasets} {

        global $this-datasets
        set $this-datasets $datasets

              set w .ui[modname]

        if [ expr [winfo exists $w.f.options] ] {

            for {set i 0} {$i < 2} {incr i 1} {
          if [ expr [winfo exists $w.f.options.datasets.$i] ] {
              pack forget $w.f.options.datasets.$i
          }
            }

            set i 0

            foreach dataset $datasets {
          if [ expr [winfo exists $w.f.options.datasets.$i] ] {
              $w.f.options.datasets.$i configure -text $dataset
          } else {
              set len [expr [string length $dataset] + 5 ]
              sci_label $w.f.options.datasets.$i -text $dataset \
            -anchor w -just left -fore darkred 
          }

          pack $w.f.options.datasets.$i -side top -anchor nw -fill x

          incr i 1
            }
        }
    }
}


