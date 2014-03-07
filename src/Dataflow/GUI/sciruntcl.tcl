puts "Load Tcl 8.6"
package require Tcl 8.6
puts "Load Tk 8.6"
package require Tk 8.6
puts "Load itcl 4.0"
package require itcl 4.0
puts "Load itk 4.0"
package require itk 4.0
puts "Load Iwidgets 4.0.1"
package require Iwidgets 4.0.1
puts "Load BLT 2.4"
package require BLT 2.4
puts "Done configuring TCL"

namespace eval ::sciruntcl {
    namespace export *

    variable library [file dirname [info script]]
    variable version 1.0
}

lappend auto_path $sciruntcl::library 
package provide sciruntcl $sciruntcl::version

source [file dirname [info script]]/defaults.tcl
source [file dirname [info script]]/Module.tcl
#source [netedit getenv SCIRUN_SRCDIR]/Dataflow/GUI/Connection.tcl
#source [netedit getenv SCIRUN_SRCDIR]/Dataflow/GUI/Port.tcl
#source [netedit getenv SCIRUN_SRCDIR]/Dataflow/GUI/Subnet.tcl
#source [netedit getenv SCIRUN_SRCDIR]/Dataflow/GUI/UIvar.tcl
#source [netedit getenv SCIRUN_SRCDIR]/Dataflow/GUI/Range.tcl
#source [netedit getenv SCIRUN_SRCDIR]/Dataflow/GUI/ToolTipText.tcl
#source [netedit getenv SCIRUN_SRCDIR]/Dataflow/GUI/NetworkEditor.tcl
