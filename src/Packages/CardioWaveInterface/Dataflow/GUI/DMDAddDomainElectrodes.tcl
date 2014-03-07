itcl::class CardioWaveInterface_DiscreteMultiDomain_DMDAddDomainElectrodes {
    inherit Module
     constructor { {args ""} } {
        eval configure $args
        set name DMDAddDomainElectrodes
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
      wm minsize $w 100 50

      iwidgets::labeledframe $w.m -labeltext "DOMAIN ELECTRODE"
      set ref [$w.m childsite]
      pack $w.m -fill both -expand yes


      radiobutton $ref.single -text "Use single domain:" \
            -variable $this-range -value "0"

      label $ref.lab1 -text "Domain"
      label $ref.lab2 -text "Min Domain"
      label $ref.lab3 -text "Max Domain"
      entry $ref.en1 -textvariable $this-electrodedomain -width 4 
            
      radiobutton $ref.range -text "Use domain range (min,max):" \
            -variable $this-range -value "1"

      entry $ref.en2 -textvariable $this-electrodedomainmin -width 4
      entry $ref.en3 -textvariable $this-electrodedomainmax -width 4
      
      grid $ref.single -row 0 -column 0 -columnspan 4 -sticky w
      grid $ref.en1 -row 1 -column 1 -sticky w
      grid $ref.lab1 -row 1 -column 0 -sticky w
      grid $ref.range -row 2 -column 0 -columnspan 4 -sticky w
      grid $ref.en2 -row 3 -column 1 -sticky w
      grid $ref.lab2 -row 3 -column 0 -sticky w
      grid $ref.en3 -row 3 -column 3 -sticky w
      grid $ref.lab3 -row 3 -column 2 -sticky w

      makeSciButtonPanel $w $w $this
      moveToCursor $w

    }
}


