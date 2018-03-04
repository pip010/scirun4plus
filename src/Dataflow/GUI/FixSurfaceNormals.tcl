#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

itcl::class SCIRun_ChangeMesh_FixSurfaceNormals {
  inherit Module
    constructor { {args ""} } {
      eval configure $args
      set name FixSurfaceNormals
      set_defaults
  }

  method set_defaults {} {
  }

  method ui {} {
		set w .ui[modname]
		if {[winfo exists $w]} {
			raise $w
			return;
		}

		sci_toplevel $w
		#wm minsize $w 80 130


		sci_frame $w.form -relief groove -borderwidth 2
		sci_label $w.form.t1 -text "Fix surface normals"
		pack $w.form.t1 

			
		sci_checkbutton $w.form.output_elems -text "Output Inverted" \
			-variable $this-output_inverted

		sci_frame $w.form.f 
		
		sci_label $w.form.f.label -text "Seed (index of element):"
		sci_entry $w.form.f.value -textvariable $this-seed
		pack $w.form.f.label $w.form.f.value -side left -anchor n
				
		pack $w.form.output_elems  $w.form.f -side top -anchor w

		pack $w.form -side top -e y -f both -padx 5 -pady 5

		makeSciButtonPanel $w $w $this
		moveToCursor $w
  }
}
