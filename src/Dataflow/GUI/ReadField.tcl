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


# GUI for ReadField module
# by Samsonov Alexei
# December 2000

catch {rename SCIRun_DataIO_ReadField ""}

itcl::class SCIRun_DataIO_ReadField {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
	set name ReadField
    }

    method ui {} {
        global $this-filename

	set w .ui[modname]

	if {[winfo exists $w]} {
	    return
	}

	sci_toplevel $w -class TkFDialog
	# place to put preferred data directory
	# it's used if $this-filename is empty

        set initdir ""

        if {[info exists env(SCIRUN_DATA)]} {
              set initdir $env(SCIRUN_DATA)
        } elseif {[info exists env(SCI_DATA)]} {
              set initdir $env(SCI_DATA)
        } elseif {[info exists env(PSE_DATA)]} {
              set initdir $env(PSE_DATA)
        } elseif {[info exists env(SCIRUN_OBJDIR)]} {
              set initdir $env(SCIRUN_OBJDIR)
        }
	
	#######################################################
	# to be modified for particular reader

	# extansion to append if no extension supplied by user
	set defext ".fld"
	set title "Open field file"
	
	######################################################
	
	# Unwrap $this-types into a list.
	set tmp1 [set $this-types]
	set tmp2 [eval "set tmp3 $tmp1"]

	makeOpenFilebox \
	    -parent $w \
	    -filevar $this-filename \
	    -filename_basevar $this-filename_base \
	    -setcmd "wm withdraw $w" \
	    -command "$this-c needexecute" \
	    -cancel "wm withdraw $w" \
	    -title $title \
	    -filetypes $tmp2 \
	    -initialdir $initdir \
	    -defaultextension $defext \
	    -allowMultipleFiles $this \
	    -selectedfiletype $this-filetype \
            -fromenv $this-from-env \
            -delayvar $this-delay
	
	moveToCursor $w	
    }
}
