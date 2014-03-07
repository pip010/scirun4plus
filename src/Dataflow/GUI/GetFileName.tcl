#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  License for the specific language governing rights and limitations under
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


# GUI for SCIRun_String_GetFileName module

# This GUI interface is for selecting a file name via the makeOpenFilebox.

itcl::class SCIRun_String_GetFileName {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
        set name GetFileName
    }

    method ui {} {

        global $this-filename
        
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }

        sci_toplevel $w -class TkFDialog

        set initdir ""
      
        # place to put preferred data directory
        # it's used if $this-filename is empty
      
        # Use the standard data dirs
        # I guess there is no .mat files in there
        # at least not yet

        if {[info exists env(SCIRUN_DATA)]} {
              set initdir $env(SCIRUN_DATA)
        } elseif {[info exists env(SCI_DATA)]} {
              set initdir $env(SCI_DATA)
        } elseif {[info exists env(PSE_DATA)]} {
              set initdir $env(PSE_DATA)
        } elseif {[info exists env(SCIRUN_OBJDIR)]} {
              set initdir $env(SCIRUN_OBJDIR)
        }
      
        makeOpenFilebox \
            -parent $w \
            -filevar $this-filename \
	     -setcmd "wm withdraw $w" \
            -command "$this-c sendfilename" \
            -commandname "Execute" \
            -cancel "wm withdraw $w" \
            -title "Select file" \
            -filetypes {{ "All files" "*.*" } }\
            -initialdir $initdir \
            -defaultextension "*.*" \
            -allowMultipleFiles $this \
            -selectedfiletype 0 \
	    -filename_basevar $this-filename_base \
	    -numbervar $this-number_in_series \
	    -delayvar $this-delay \
	    -pinnedvar $this-pinned
	
        moveToCursor $w
    }
}
