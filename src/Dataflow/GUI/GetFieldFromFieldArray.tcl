itcl::class SCIRun_FieldArray_GetFieldFromFieldArray {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name GetFieldFromFieldArray
        set_defaults
    }

    method set_defaults {} {
    }

    method update_range { } {
    
      set w .ui[modname]
      if {[winfo exists $w]} {
      
        upvar \#0 $this-index-min min $this-index-max max         
        $w.selectfield configure -from $min -to $max -variable $this-index -command "$this-c needexecute"
      }
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            return
        }
        toplevel $w

        sci_scale $w.selectfield -from 0 -to 1 -variable $this-index -showvalue true -orient horizontal -relief groove -label "Field Index:" 
        pack $w.selectfield -side top -anchor n -fill x -padx 8 -pady 5
        update_range

        makeSciButtonPanel $w $w $this
        moveToCursor $w
    }
}


