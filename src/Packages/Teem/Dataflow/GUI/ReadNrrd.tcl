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


# GUI for ReadNrrd module
# by Samsonov Alexei
# December 2000

catch {rename Teem_DataIO_ReadNrrd ""}

itcl::class Teem_DataIO_ReadNrrd {
    inherit Module

     constructor { {args ""} } {
        eval configure $args
      set name ReadNrrd
    }

    method ui {} {
      #global env
      #global $this-filename
      
      set w .ui[modname]
      
      if {[winfo exists $w]} {
          return
      }
      
      sci_toplevel $w -class TkFDialog
      
      # place to put preferred data directory
      # it's used if $this-filename is empty
      set initdir [netedit getenv SCIRUN_DATA]
      
      #######################################################
      # to be modified for particular reader
      
      # extansion to append if no extension supplied by user
      set defext ".nrrd"
      set title "Open nrrd file"

      set tmp1 [set $this-types]
      set tmp2 [eval "set tmp3 $tmp1"]
      
      # file types to appers in filter box
      #set types {
      #      {{Nrrd Files}         {.nhdr .nrrd .png .txt .vtk}   }
	    #{{NrrdData File}      {.nd}           }
	    #{{VFF File}           {.vff}          }
	    #{{PICT File}          {.pic .pict}    }
	    #{{Vista File}         {.v}    }
	    #{{All Files}          {.*}            }
      #}
      
      ######################################################
	
      makeOpenFilebox \
          -parent $w \
          -filevar $this-filename \
          -setcmd "wm withdraw $w" \
          -command "$this-c needexecute" \
          -cancel "wm withdraw $w" \
          -title $title \
          -filetypes $tmp2 \
          -initialdir $initdir \
          -defaultextension $defext \
          -selectedfiletype $this-filetype \
          -fromenv $this-from-env

      moveToCursor $w
    }
}
