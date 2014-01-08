itcl::class CardioWaveInterface_DiscreteMultiDomain_DMDAddStimulusSeries {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name DMDAddStimulusSeries
        set_defaults
    }

    method set_defaults {} {
    }

    method ui {} {
      
      set w .ui[modname]
      if {[winfo exists $w]} {
          return
      }

      toplevel $w 
      wm minsize $w 100 150

      iwidgets::labeledframe $w.m -labeltext "BLOCK PULSE STIMULUS"
      set stim [$w.m childsite]
      pack $w.m -fill both -expand yes

      radiobutton $stim.lab1 -text "Stimulus Domain Single" -variable $this-stim-domain-range -value "0"
      radiobutton $stim.lab1a -text "Stimulus Domain Range" -variable $this-stim-domain-range -value "1"

      label $stim.lab2 -text "Stimulus Current (mA)"
      label $stim.lab3 -text "Stimulus Start (ms)"
      label $stim.lab4 -text "Stimulus End (ms)"
      label $stim.lab5 -text "Stimulus Repeat (ms)"

      entry $stim.en1 -textvariable $this-stim-domain -width 8
      entry $stim.en1a -textvariable $this-stim-domain-min -width 8
      entry $stim.en1b -textvariable $this-stim-domain-max -width 8
      entry $stim.en2 -textvariable $this-stim-current -width 8
      entry $stim.en3 -textvariable $this-stim-start -width 8
      entry $stim.en4 -textvariable $this-stim-end -width 8
      entry $stim.en5 -textvariable $this-stim-duration -width 8
	
      checkbutton $stim.cb1 -text "Current is Current Density" -variable $this-stim-is-current-density
      checkbutton $stim.cb2 -text "Define stimulus only for elements contained within geometry" -variable $this-stim-useelements
      grid $stim.lab1 -row 0 -column 0 -sticky w
      grid $stim.lab1a -row 1 -column 0 -sticky w
      grid $stim.lab2 -row 2 -column 0 -sticky w
      grid $stim.lab3 -row 3 -column 0 -sticky w
      grid $stim.lab4 -row 4 -column 0 -sticky w
      grid $stim.lab5 -row 5 -column 0 -sticky w
      grid $stim.en1 -row 0 -column 1 -sticky w
      grid $stim.en1a -row 1 -column 1 -sticky w
      grid $stim.en1b -row 1 -column 2 -sticky w
      grid $stim.en2 -row 2 -column 1 -sticky w
      grid $stim.en3 -row 3 -column 1 -sticky w
      grid $stim.en4 -row 4 -column 1 -sticky w
      grid $stim.en5 -row 5 -column 1 -sticky w
      grid $stim.cb1 -row 6 -column 0 -columnspan 3 -sticky w      
      grid $stim.cb2 -row 7 -column 0 -columnspan 3 -sticky w     


      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }

}


