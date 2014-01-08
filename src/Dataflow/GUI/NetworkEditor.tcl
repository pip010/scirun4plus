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

set SCIRUN_SRCDIR [netedit getenv SCIRUN_SRCDIR]
set smallIcon [image create photo -file "$SCIRUN_SRCDIR/pixmaps/scirun-icon-small.ppm"]
set splashImageFile "$SCIRUN_SRCDIR/Applications/SCIRun/scisplash.ppm"
set splashAboutFile "$SCIRUN_SRCDIR/Applications/SCIRun/sciabout.ppm"
set bioTensorSplashImageFile "$SCIRUN_SRCDIR/Packages/Teem/Dataflow/GUI/splash-tensor.ppm"
set bioFEMSplashImageFile "$SCIRUN_SRCDIR/Packages/BioPSE/Dataflow/GUI/splash-biofem.ppm"
set bioImageSplashImageFile "$SCIRUN_SRCDIR/Packages/Teem/Dataflow/GUI/splash-bioimage.ppm"
set fusionViewerSplashImageFile "$SCIRUN_SRCDIR/Packages/Fusion/Dataflow/GUI/splash-fusionviewer.ppm"
set levelSetSegmenterViewerSplashImageFile "$SCIRUN_SRCDIR/Applications/SCIRun/scisplash.ppm"
set padlockImageFile "$SCIRUN_SRCDIR/main/padlock.ppm"

set modname_font "-Adobe-Helvetica-Medium-R-Normal-*-12-120-75-*"
set ui_font "-Adobe-Helvetica-Medium-R-Normal-*-12-120-75-*"
set time_font "-Adobe-Courier-Medium-R-Normal-*-12-120-75-*"

set small_font_ui "-Adobe-Courier-Medium-R-Normal-*-12-120-75-*"
set small_font_mod "-Adobe-Helvetica-Medium-R-Normal-*-12-120-75-*"

set medium_font_ui "-Adobe-Courier-Bold-R-Normal-*-14-120-75-*"
set medium_font_mod "-Adobe-Helvetica-Bold-R-Normal-*-14-120-75-*"

set firstIcon 1

set mainCanvasWidth    4500.0
set mainCanvasHeight   4500.0
set maincanvas .frame.editor.canvas
set minicanvas .frame.editor.canvas.mini

# Records mouse position at button press to bring up menus at 
set mouseX 0
set mouseY 0

set Subnet(Subnet0_minicanvas) $minicanvas
set Subnet(Subnet0_canvas) $maincanvas
set Subnet(Subnet0_Name) Main
set Subnet(Subnet0_Modules) ""
set Subnet(Subnet0_connections) ""
set Subnet(Subnet0_netversion) 0
set Subnet(Subnet0) 0
set Subnet(Loading) 0
set Subnet(num) 0

set inserting 0
set netedit_savefile ""
set NetworkChanged 0
set CurrentlySelectedModules ""
set disable_network_locking "0"
set network_executing "0"
trace variable network_executing w handle_network_executing

# Variables used to update the progress "loading SCIRun" display.
set totalProgressSteps  0
set currentProgressStep 0
set currentPercent      0; # Integer range 0-100

set ENV_SCIRUN_DATA [netedit getenv SCIRUN_DATA]
set ENV_SCIRUN_MYDATA_DIR [netedit getenv SCIRUN_MYDATA_DIR]
set ENV_SCIRUN_NET_DIR [netedit getenv SCIRUN_NET_DIR]
set ENV_SCIRUN_DATASET [netedit getenv SCIRUN_DATASET]

set ENV_SCIRUN_NET_SUBSTITUTE_DATADIR [netedit getenv SCIRUN_NET_SUBSTITUTE_DATADIR]
set ENV_SCIRUN_NET_RELATIVE_FILENAMES [netedit getenv SCIRUN_NET_RELATIVE_FILENAMES]

set ENV_SCIRUN_LARGE_MODULES [netedit getenv SCIRUN_LARGE_MODULES]
set ENV_SCIRUN_STRAIGHT_COONECTIONS [netedit getenv SCIRUN_STRAIGHT_CONNECTIONS]

set ENV_SCIRUN_STARTMATLAB [netedit getenv SCIRUN_STARTMATLAB]
set ENV_SCIRUN_MATLABTIMEOUT [netedit getenv SCIRUN_MATLABTIMEOUT]

proc set_network_executing {val} {
    global network_executing disable_network_locking
    if { $disable_network_locking == "1" } {
	set network_executing "0"
	return
    }
    set network_executing $val
}

proc restore_not_executing_interface {} {
}

proc disable_netedit_locking {} {
}

proc handle_network_executing { var op1 op2} {
}

proc setIcons { { w . } { size small } } {
    global firstIcon tcl_platform
    set srcdir [netedit getenv SCIRUN_SRCDIR]
    if { [string length $srcdir] } {
	set bitmap $srcdir/pixmaps/scirun-icon-$size.xbm
	set inverted $srcdir/pixmaps/scirun-icon-$size-inverted.xbm
    
        global OnWindows
        if { $OnWindows && $size == "small" && $firstIcon == 1 } {
            # make all subsequent windows use this icon.  This prevents an annoying problem
            # where an empty window pops up and then fills up and moves to the proper location
            wm iconbitmap $w -default @$inverted
            set firstIcon 0
	} elseif { ! $OnWindows } {
            wm iconbitmap $w @$inverted	
	}
	wm iconmask $w @$bitmap
    }
}


# envBool <variable name>
#
#   envBool will query the enviromnent varaible's value and return as a boolean
#   Usage example:  envBool SCIRUN_INSERT_NET_COPYRIGHT
#
#   Turns 'true/false', 'on/off', 'yes/no', '1/0' into '1/0' respectively
#   (Non-existent and empty  variables are treated as 'false'.)
#
#   This function is case insensitive.
#
proc envBool { var } {
    set val [netedit getenv $var] ; # returns blank string if variable not set
    if { [string equal $val ""] } {
      return 0; # blank value is taken to mean false
    }
    if { ![string is boolean $val] } {
      puts "TCL envBool: Cannot determine boolean value of env: $var=$val"
      return 1; # follows the C convention of any non-zero value equals true
    }
    return [string is true $val]
}

proc safeSetWindowGeometry { w geom } {
    set realgeom [split [wm geometry $w] +]
    set geom [split $geom +]
    set sizepos [lsearch $geom *x*]
    if { $sizepos == -1 } {
      set width [winfo width $w]
      set height [winfo height $w]
    } else {
      set size [split [lindex $geom $sizepos] x]
      set geom [lreplace $geom $sizepos $sizepos]
      set width [lindex $size 0]
      set height [lindex $size 1]
    }

    set xoff [lindex $realgeom 1]
    set yoff [lindex $realgeom 2]

    if { [llength $geom] } {
      set xoff [lindex $geom end-1]
      set yoff [lindex $geom end]
    }
    
    set xoff    [expr ($xoff < 0) ? 0 : $xoff]
    set yoff    [expr ($yoff < 0) ? 0 : $yoff]

    wm geometry $w ${width}x${height}+${xoff}+${yoff}
}

rename toplevel __TCL_toplevel__
proc toplevel { args } {
    set win [uplevel 1 __TCL_toplevel__ $args]
    setIcons [lindex $args 0]
    set varname [string range $win 3 end]-ui_geometry
    upvar \#0 $varname geometry
    if { [info exists geometry] } {
	safeSetWindowGeometry $win $geometry
    }
    return $win
}


proc geometryTrace { args } {
    global geometry Subnet
    if { [info exists geometry] } {
	set w .
	if { $Subnet(Loading) } {
	    set w .subnet$Subnet(Loading)
	}
	safeSetWindowGeometry $w $geometry
    }
}

proc makeNetworkEditor {} {

    global Color maincanvas minicanvas mainCanvasHeight mainCanvasWidth
    global geometry OnWindows
    global ui_font modname_font time_font

    set ostype [netedit getenv OS]
    if { $ostype == "Windows_NT" } {
        set OnWindows 1
    } else {
        set OnWindows 0
    }

    wm protocol . WM_DELETE_WINDOW { NiceQuit }
    wm minsize . 150 150

    trace variable geometry w geometryTrace

    set geom [netedit getenv SCIRUN_GEOMETRY]
    if { [string length $geom] } {
      safeSetWindowGeometry . $geom
    } else {
      safeSetWindowGeometry . 800x800
    }

    wm title . "SCIRun"
    setIcons . large

    loadToolTipText

# BUILD THE TOP MENU BAR

    frame .main_menu -borderwidth 0 -background $Color(MenuBarBackGround)
    
    menubutton .main_menu.file -text "File" -underline 0 \
      -menu .main_menu.file.menu -bd 0\
      -background $Color(MenuBarBackGround) \
      -foreground $Color(MenuBarForeGround) \
      -activebackground $Color(MenuSelectBackGround)  \
      -activeforeground $Color(MenuSelectForeGround)

    menu .main_menu.file.menu -tearoff false -bd 2 -activeborderwidth 0 \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(MenuBackGround) \
     -foreground $Color(MenuForeGround)

    # Create the "File" Menu sub-menus.  Create (most of) them in the
    # disabled state.  They will be enabled when all packages are loaded.
    .main_menu.file.menu add command -label "Load..." -underline 0 \
      -command "popupLoadMenu" -state disabled
    .main_menu.file.menu add command -label "Insert..." -underline 0 \
      -command "popupInsertMenu" -state disabled
    .main_menu.file.menu add command -label "Save" -underline 0 \
      -command "popupSaveMenu" -state disabled
    .main_menu.file.menu add command -label "Save As..." -underline 0 \
      -command "popupSaveAsMenu" -state disabled

    .main_menu.file.menu add separator

    .main_menu.file.menu add command -label "Clear Network" -underline 0 \
      -command "ClearCanvas" -state disabled
    .main_menu.file.menu add command -label "Select All" -underline 0 \
      -command "selectAll" -state disabled
    .main_menu.file.menu add command -label "Execute All" -underline 0 \
      -command "updateRunDateAndTime 0; netedit scheduleall" -state disabled

    .main_menu.file.menu add separator
      .main_menu.file.menu add command -label "Create Module Skeleton..." \
        -underline 0 -command "ComponentWizard" -state disabled
      .main_menu.file.menu add separator
    
# This was added by Mohamed Dekhil to add some infor to the net
#    .main_menu.file.menu add command -label "Network Properties" -underline 0 \
#      -command "popupInfoMenu"

#    .main_menu.file.menu add separator
#    .main_menu.file.menu add command -label "Edit Configuration" -underline 0 \
#      -command "editrcfile"

#    .main_menu.file.menu add separator
    .main_menu.file.menu add command -label "Quit" -underline 0 \
	    -command "NiceQuit"

    menubutton .main_menu.help -text "Help" -underline 0 \
      -menu .main_menu.help.menu -direction below -bd 0 \
      -background $Color(MenuBarBackGround) \
      -foreground $Color(MenuBarForeGround) \
      -activebackground $Color(MenuSelectBackGround)  \
      -activeforeground $Color(MenuSelectForeGround)
      
    menu .main_menu.help.menu -tearoff false -bd 2 -activeborderwidth 0 \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(MenuBackGround) \
     -foreground $Color(MenuForeGround)
    .main_menu.help.menu add checkbutton -label "Show Tooltips" -underline 0 \
      -variable tooltipsOn

    # Mac hack to fix size of 'About' window ... sigh... 
    .main_menu.help.menu add command -label "About..." -underline 0 \
      -state disabled -command  "showProgress 2 none 1"

    .main_menu.help.menu add command -label "Credits..." -underline 0 \
      -state disabled -command  "showProgress 1 none 1"

    .main_menu.help.menu add command -label "License..." -underline 0 \
      -command  "licenseDialog" -state disabled

    global ModuleFilter
    set ModuleFilter ""
    label .main_menu.highlightlabel -text "Show Module:" -background $Color(MenuBarBackGround)
    sci_entry .main_menu.highlight \
      -textvariable ModuleFilter -width 26

    label .main_menu.disclaimer \
      -text "SCIRun v[netedit getenv SCIRUN_VERSION] (revision [netedit getenv SCIRUN_SVN_REVISION])" \
      -background $Color(MenuBarBackGround) \
      -foreground $Color(SCIText)

    menubutton .main_menu.module -text "Modules" -underline 0 \
      -menu .main_menu.module.menu -bd 0\
      -background $Color(MenuBarBackGround) \
      -foreground $Color(MenuBarForeGround) \
      -activebackground $Color(MenuSelectBackGround)  \
      -activeforeground $Color(MenuSelectForeGround)
     
    menu .main_menu.module.menu -tearoff false -bd 2 -activeborderwidth 0 \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(MenuBackGround) \
     -foreground $Color(MenuForeGround) \
     -postcommand UpdateHighLightsPackageMenu

    menubutton .main_menu.subnet -text "Subnets" -underline 0 \
      -menu .main_menu.subnet.menu -bd 0 \
      -background $Color(MenuBarBackGround) \
      -foreground $Color(MenuBarForeGround) \
      -activebackground $Color(MenuSelectBackGround)  \
      -activeforeground $Color(MenuSelectForeGround)

    menu .main_menu.subnet.menu -tearoff false -postcommand createSubnetMenu \
     -bd 2 -activeborderwidth 0 \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(MenuBackGround) \
     -foreground $Color(MenuForeGround)

    menubutton .main_menu.toolkits -text "Toolkits" -underline 0 \
      -menu .main_menu.toolkits.menu -bd 0 \
      -background $Color(MenuBarBackGround) \
      -foreground $Color(MenuBarForeGround) \
      -activebackground $Color(MenuSelectBackGround)  \
      -activeforeground $Color(MenuSelectForeGround)

    menu .main_menu.toolkits.menu -tearoff false -postcommand createToolkitsMenu \
     -bd 2 -activeborderwidth 0 \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(MenuBackGround) \
     -foreground $Color(MenuForeGround)

    pack .main_menu.file -side left
    pack .main_menu.module -side left
    pack .main_menu.subnet -side left
    pack .main_menu.toolkits -side left
    pack .main_menu.help -side left

    pack .main_menu.disclaimer -side right -padx 6
    pack .main_menu.highlight -side right 
    pack .main_menu.highlightlabel -side right
    
#    tk_menuBar .main_menu .main_menu.file

# BUILD SCIRUN MENUBAR
    
    frame .scirun_menu -background $Color(MainBackGround) 
    
    sci_button .scirun_menu.execute_all \
      -command "updateRunDateAndTime 0; netedit scheduleall" -text "Execute All" 
    
    frame .scirun_menu.space -width 24 -background $Color(MainBackGround)
    pack .scirun_menu.execute_all  -padx 8 -side left

    label .scirun_menu.indicator -text "0/0" -background $Color(MainBackGround) -foreground black
    sci_button .scirun_menu.configure -text "< Configure >" -width 12
    sci_button .scirun_menu.detach -text "+"  

    bind .scirun_menu.configure <ButtonPress-1> "SwitchConfigMenu3"
    bind .scirun_menu.configure <ButtonPress-2> "SwitchConfigMenu"
    bind .scirun_menu.configure <ButtonPress-3> "SwitchConfigMenu2"
    bind .scirun_menu.detach <ButtonPress-1>    "SwitchConfigMenu2"
    
    pack .scirun_menu.space .scirun_menu.detach .scirun_menu.configure .scirun_menu.indicator -side right -padx 1 -pady 3
    
    frame .scirun_menu.inset -relief sunken -height 16 \
                -borderwidth 1 -width [expr 80-4] -background $Color(Trough)
    pack .scirun_menu.inset -side right -fill y -padx 1 -pady 3
    frame .scirun_menu.inset.graph -relief raised \
                -width 0 -borderwidth 1 -background $Color(ModuleProgress)

# BUILD MAIN CANVAS

    frame .frame 
    frame .frame.editor -relief sunken -borderwidth 2 -bg $Color(Basecolor) 

    canvas $maincanvas -bg "$Color(NetworkEditor)" \
        -scrollregion "0 0 $mainCanvasWidth $mainCanvasHeight"
    pack $maincanvas -expand 1 -fill both
  
    # bgRect is just a rectangle drawn on the neteditFrame Canvas
    # so that the Modules List Menu can be bound to it using mouse
    # button 3.  The Modules List Menu can't be bound to the canvas
    # itself because mouse events are sent to both the objects on the
    # canvas (such as the lines connection the modules) and the canvas.
     eval $maincanvas create rectangle [$maincanvas cget -scrollregion] \
      -fill "$Color(NetworkEditor)" -outline "$Color(NetworkEditor)" \
      -tags bgRect 

    $maincanvas configure \
      -xscrollcommand "updateViewAreaBox; .frame.hscroll set" \
      -yscrollcommand "updateViewAreaBox; .frame.vscroll set"

    scrollbar .frame.hscroll -orient horizontal \
      -command "$maincanvas xview" -elementborderwidth 1 \
      -bd 1 -troughcolor $Color(Trough) \
      -background $Color(MainBackGround) 
    
    scrollbar .frame.vscroll \
      -command "$maincanvas yview" -elementborderwidth 1 \
      -bd 1 -troughcolor $Color(Trough) \
      -background $Color(MainBackGround) 

    # Layout the scrollbars and canvas in the bottom pane
    grid .frame.editor .frame.vscroll .frame.hscroll
    grid columnconfigure .frame 0 -weight 1 
    grid rowconfigure    .frame 0 -weight 1 
    grid config .frame.editor -column 0 -row 0 \
	    -columnspan 1 -rowspan 1 -sticky "snew" 
    grid config .frame.hscroll -column 0 -row 1 \
	    -columnspan 1 -rowspan 1 -sticky "ew" -pady 2
    grid config .frame.vscroll -column 1 -row 0 \
	    -columnspan 1 -rowspan 1 -sticky "sn" -padx 0

    pack .main_menu -fill x
    pack .scirun_menu -fill x
    pack .frame -expand 1 -fill both -after .main_menu

    # Create Mini Network Editor
    global miniCanvasWidth miniCanvasHeight
    canvas $minicanvas -bg $Color(NetworkEditorSmall) -width 120 -height 120 -borderwidth 1 -relief raised
    pack $minicanvas -anchor se -side bottom -padx 3 -pady 3

    $minicanvas create rectangle 0 0 1 1 -outline white -width 2 -tag "viewAreaBox"
    initInfo
    
    # SETUP THE SETTING WINDOWS CONFIGURATIONS

    # DEFINE THE WINDOW AND AN EMPTY FRAME INSIDE
    toplevel .detached
    frame .detached.f -background $Color(MainBackGround)
    pack .detached.f -side top -anchor n -fill y -expand yes

    wm title .detached "SCIRun Configuration"
    wm sizefrom .detached user
    wm positionfrom .detached user
    wm protocol .detached WM_DELETE_WINDOW "RemoveConfigFrame"
    wm withdraw .detached
    
    # This is the frame for the geometry controls
    iwidgets::scrolledframe .configframe -width 350 -height 230 \
      -hscrollmode dynamic\
      -vscrollmode none\
      -background $Color(MainBackGround) \
      -foreground $Color(MainForeGround) -borderwidth 0
      
    set sc [.configframe component horizsb]  
    $sc configure -elementborderwidth 1 \
      -bd 1 -troughcolor $Color(Trough) \
      -background $Color(MainBackGround) 

    set sc [.configframe component vertsb] 
    $sc configure -elementborderwidth 1 \
      -bd 1 -troughcolor $Color(Trough) \
      -background $Color(MainBackGround) 
        
    # get the childsite to add stuff to
    set configframe [.configframe childsite]

    frame $configframe.f
    pack $configframe.f -side top -anchor nw -fill x -expand yes
    pack $configframe -side top -anchor nw -padx 10 -fill x -expand yes

    global IsAttached IsDisplayed

    set IsAttached 1
    set IsDisplayed 0
    
    InitConfigFrame .detached.f $configframe.f
    InitConfigFrame $configframe.f .detached.f
    
#    wm withdraw .

}

proc RemoveConfigFrame {} {
    global IsAttached IsDisplayed

    if { $IsAttached != 0 } {
        pack forget .configframe
        set height [expr [winfo height .]-[winfo height .configframe]]
        wm geometry . [winfo width .]x${height}
        update
    } else { 
        wm withdraw .detached
    }
    
    .scirun_menu.configure configure -text "< Configure >"
    .scirun_menu.detach configure -text "+"
    bind .scirun_menu.configure <ButtonPress-1> "SwitchConfigMenu3"
    bind .scirun_menu.configure <ButtonPress-2> "SwitchConfigMenu"
    bind .scirun_menu.configure <ButtonPress-3> "SwitchConfigMenu2"
    bind .scirun_menu.detach <ButtonPress-1>    "SwitchConfigMenu2"
   
    set IsDisplayed 0
  }

proc AddConfigFrame {} {
    global IsAttached IsDisplayed
    
    if { $IsAttached!=0} {
        pack .configframe -anchor w -side bottom -before .scirun_menu -fill x
        set w1 [winfo width .]
        set w2 [winfo width .configframe]
        set width [expr $w1 > $w2 ? $w1 : $w2]
        set height [expr [winfo height .]+[winfo reqheight .configframe]]
        wm geometry . ${width}x${height}
        update
    } else {
        wm deiconify $detachedFr
    }
    
    .scirun_menu.configure configure -text "> Configure <"
    .scirun_menu.detach configure -text "S"
    bind .scirun_menu.configure <ButtonPress-1> "RemoveConfigFrame"
    bind .scirun_menu.configure <ButtonPress-2> "SwitchConfigMenu"
    bind .scirun_menu.configure <ButtonPress-3> "RemoveConfigFrame"
    bind .scirun_menu.detach <ButtonPress-1>     "SwitchConfigMenu2"
   
    set IsDisplayed 1
  }


proc InitConfigFrame {m m2} {
    if { ![winfo exists $m] } return

    global Color
    frame $m.config -background $Color(UIBackGround) -bd 1 -relief sunken
    pack $m.config -fill x -expand yes -side top

    iwidgets::tabnotebook $m.config.tabs -height 210 -width 700 -tabpos n \
      -background $Color(UIBackGround) \
      -foreground $Color(UIForeGround) \
      -backdrop $Color(UIBackDrop) \
      -tabbackground $Color(UIBackDrop) \
      -bevelamount 0 -gap 0 -borderwidth 0 \
      -equaltabs false
      
    $m.config.tabs add -label "Network Editor"
    $m.config.tabs add -label "Options"
    $m.config.tabs add -label "Log"
    $m.config.tabs select 0

    pack $m.config.tabs -fill x -expand yes -anchor w

    set general   [$m.config.tabs childsite 1]
    set module    [$m.config.tabs childsite 0]
    set log       [$m.config.tabs childsite 2]

    sci_frame $module.f1
    pack $module.f1 -side left -anchor nw
    
    sci_labeledframe $module.f1.f -labeltext "Layout Options"
    pack $module.f1.f -side top -anchor nw -pady 5 -padx 5
    set f [$module.f1.f childsite]    
    
    global ENV_SCIRUN_LARGE_MODULES
    sci_checkbutton $f.lmodules -text "Use large modules in editor."\
      -variable ENV_SCIRUN_LARGE_MODULES \
      -command "global ENV_SCIRUN_LARGE_MODULES; netedit update_env SCIRUN_LARGE_MODULES \$ENV_SCIRUN_LARGE_MODULES; UpdateCanvas"
    grid $f.lmodules -row 0 -column 0 -sticky w -pady 5 -padx 5 

    global ENV_SCIRUN_STRAIGHT_CONNECTIONS
    sci_checkbutton $f.straightlines -text "Use straight dataflow pipelines in editor."\
      -variable ENV_SCIRUN_STRAIGHT_CONNECTIONS \
      -command "global ENV_SCIRUN_STRAIGHT_CONNECTIONS; netedit update_env SCIRUN_STRAIGHT_CONNECTIONS \$ENV_SCIRUN_STRAIGHT_CONNECTIONS; UpdateCanvas"
    grid $f.straightlines -row 1 -column 0 -sticky w -pady 5 -padx 5 

    sci_frame $general.f1
    sci_frame $general.f2
    pack $general.f1 $general.f2 -side left -anchor nw

    sci_labeledframe $general.f1.f -labeltext "SCIRun Path Settings"
    pack $general.f1.f -side top -anchor nw -pady 5 -padx 5
    set f [$general.f1.f childsite]
    
    global ENV_SCIRUN_DATA
    sci_label $f.lab1 -text "SCIRun Data"
    sci_entry $f.ent1 -textvariable ENV_SCIRUN_DATA -width 25
    sci_button $f.but1 -text "Set" -command "global ENV_SCIRUN_DATA; netedit update_env SCIRUN_DATA \$ENV_SCIRUN_DATA"

    global ENV_MYDATA_DIR
    sci_label $f.lab2 -text "User Data"
    sci_entry $f.ent2 -textvariable ENV_MYDATA_DIR -width 25
    sci_button $f.but2 -text "Set" -command "global ENV_MYDATA_DIR;netedit update_env MYDATA_DIR \$ENV_MYDATA_DIR"
    global ENV_SCIRUN_NET_DIR
    sci_label $f.lab3 -text "SCIRun Nets"
    sci_entry $f.ent3 -textvariable ENV_SCIRUN_NET_DIR -width 25
    sci_button $f.but3 -text "Set" -command "global ENV_SCIRUN_NET_DIR;netedit update_env SCIRUN_NET_DIR \$ENV_SCIRUN_NET_DIR"
    
    grid $f.lab1 -row 0 -column 0 -sticky w -pady 1 -padx 5
    grid $f.ent1 -row 0 -column 1 -pady 1 -padx 5
    grid $f.but1 -row 0 -column 2 -pady 1 -padx 5

    grid $f.lab2 -row 2 -column 0 -sticky w -pady 1 -padx 5
    grid $f.ent2 -row 2 -column 1 -pady 1 -padx 5
    grid $f.but2 -row 2 -column 2 -pady 1 -padx 5

    grid $f.lab3 -row 1 -column 0 -sticky w -pady 1 -padx 5
    grid $f.ent3 -row 1 -column 1 -pady 1 -padx 5
    grid $f.but3 -row 1 -column 2 -pady 1 -padx 5
     
    sci_labeledframe $general.f1.g -labeltext "SCIRun Data Set"
    pack $general.f1.g -side top -anchor nw -pady 5 -padx 5 -fill x
    set g [$general.f1.g childsite]
    
    global ENV_SCIRUN_DATASET
    sci_label $g.lab1 -text "Data Set"
    sci_entry $g.ent1 -textvariable ENV_SCIRUN_DATASET -width 25
    sci_button $g.but1 -text "Set" -command "global ENV_SCIRUN_DATASET; netedit update_env SCIRUN_DATASET \$ENV_SCIRUN_DATASET"
 
    grid $g.lab1 -row 0 -column 0 -sticky w -pady 5 -padx 5
    grid $g.ent1 -row 0 -column 1 -pady 1 -padx 5
    grid $g.but1 -row 0 -column 2 -pady 1 -padx 5
 
    sci_labeledframe $general.f2.f -labeltext "SCIRun Options"
    pack $general.f2.f -side top -anchor nw -pady 5 -padx 5
    set f [$general.f2.f childsite]
 
    global ENV_SCIRUN_NET_SUBSTITUTE_DATADIR
    sci_checkbutton $f.opt1 \
      -text "The names of data files in the SCIRun Data\n directory are saved relative to this one." \
      -onvalue 1 -offvalue 0 -variable ENV_SCIRUN_NET_SUBSTITUTE_DATADIR \
      -command "global ENV_SCIRUN_NET_SUBSTITUTE_DATADIR; netedit update_env SCIRUN_NET_SUBSTITUTE_DATADIR \$ENV_SCIRUN_NET_SUBSTITUTE_DATADIR"

    global ENV_SCIRUN_NET_RELATIVE_FILENAMES
    sci_checkbutton $f.opt2 \
      -text "The names of data files are saved relative\n to the path of the network file." \
      -onvalue 1 -offvalue 0 -variable ENV_SCIRUN_NET_RELATIVE_FILENAMES \
      -command "global ENV_SCIRUN_NET_RELATIVE_FILENAMES; netedit update_env SCIRUN_NET_RELATIVE_FILENAMES \$ENV_SCIRUN_NET_RELATIVE_FILENAMES"

    grid $f.opt1 -row 0 -column 0 -sticky w -pady 5 -padx 5
    grid $f.opt2 -row 1 -column 0 -sticky w -pady 5 -padx 5


    # LOG SCREEN
    sci_scrolledtext $log.text 
    pack $log.text -fill both -expand yes
    set stext [$log.text component text]

    $stext tag configure red -foreground red
    $stext tag configure blue -foreground blue
    $stext tag configure yellow -foreground "\#808000"
    $stext tag configure black -foreground "black"

}

proc AddToLog { logtext { status "info"} } {

    set status_tag "blue"
    if { "$status" == "error" } {
       set status_tag "red" 
    } elseif { "$status" == "warning" } {
       set status_tag "yellow" 
    } elseif { "$status" == "remark" } {
       set status_tag "blue" 
    }
    
    set configframe [.configframe childsite]
    set log [$configframe.f.config.tabs childsite 2]
    set text [$log.text component text] 
    $text insert end $logtext $status_tag
    $log.text yview moveto 1

    set log [.detached.f.config.tabs childsite 2]
    set text [$log.text component text] 
    $text insert end $logtext $status_tag
    $log.text yview moveto 1
}


proc SwitchConfigMenu {} {
    global Color
    global IsAttached IsDisplayed

    if { $IsDisplayed } {
      if { $IsAttached!=0} {        
        pack forget .configframe
        set hei [expr [winfo height .]-[winfo reqheight .configframe]]
        append geom [winfo width .]x${hei}
        wm geometry . $geom
        wm deiconify .detached
        set IsAttached 0
      } else {
        wm withdraw .detached
        pack .configframe -anchor w -side bottom \
            -before .scirun_menu -fill x
        set hei [expr [winfo height .]+[winfo reqheight .configframe]]
        append geom [winfo width .]x${hei}
        wm geometry . $geom
        set IsAttached 1
      }
      
      update
      .scirun_menu.configure configure -text "> Configure <"
      .scirun_menu.detach configure -text "S"
      bind .scirun_menu.configure <ButtonPress-1> "RemoveConfigFrame"
      bind .scirun_menu.configure <ButtonPress-2> "SwitchConfigMenu"
      bind .scirun_menu.configure <ButtonPress-3> "RemoveConfigFrame"    
      bind .scirun_menu.detach <ButtonPress-1> "SwitchConfigMenu"

    }
  }

proc SwitchConfigMenu2 {} {
    global Color
    global IsAttached IsDisplayed    
    
    if { !$IsDisplayed } {
      pack forget .configframe
      wm deiconify .detached
      set IsAttached 0
      set IsDisplayed 1
      update
    
      .scirun_menu.detach configure -text "S"
      .scirun_menu.configure configure -text "> Configure <"
      bind .scirun_menu.configure <ButtonPress-1> "RemoveConfigFrame"
      bind .scirun_menu.configure <ButtonPress-2> "SwitchConfigMenu"
      bind .scirun_menu.configure <ButtonPress-3> "RemoveConfigFrame"   
      bind .scirun_menu.detach <ButtonPress-1> "SwitchConfigMenu"
    } 
  }


proc SwitchConfigMenu3 {} {
    global Color
    global IsAttached IsDisplayed 

    if { !$IsDisplayed } {
      wm withdraw .detached
      pack .configframe -anchor w -side bottom \
          -before .scirun_menu -fill x
      set hei [expr [winfo height .]+[winfo reqheight .configframe]]
      append geom [winfo width .]x${hei}
      wm geometry . $geom
      set IsAttached 1
      set IsDisplayed 1
      update
    
      .scirun_menu.configure configure -text "> Configure <"
      .scirun_menu.detach configure -text "S"
      bind .scirun_menu.configure <ButtonPress-1> "RemoveConfigFrame"
      bind .scirun_menu.configure <ButtonPress-2> "SwitchConfigMenu"
      bind .scirun_menu.configure <ButtonPress-3> "RemoveConfigFrame"   
      bind .scirun_menu.detach <ButtonPress-1>    "SwitchConfigMenu"
    } 
  }



proc UpdateNetworkProgress {cnt size} {
    if { $size  > 0 } {
      set progress [expr $cnt/double($size)]
    } else {
      set progress 1
    }
    
    set width [expr int($progress*(74))]
    set graph .scirun_menu.inset.graph
    
    if {$width == 0} { 
      place forget $graph
    } else {
      $graph configure -width $width
      place $graph -relheight 1 -anchor nw
    }
    
    set progtext "$cnt/$size"
    .scirun_menu.indicator configure -text $progtext
}

proc canvasScroll { canvas { dx 0.0 } { dy 0.0 } {multiplier 1.0} } {
    if {$dx!=0.0} {$canvas xview moveto [expr $dx*$multiplier +[lindex [$canvas xview] 0]]}
    if {$dy!=0.0} {
    $canvas yview moveto [expr $dy*$multiplier +[lindex [$canvas yview] 0]]}

    raise .frame.editor.canvas.mini
}

# Activate the "File" menu items - called from C after all packages are loaded
proc activate_file_submenus { } {
    global maincanvas minicanvas    
    
    .main_menu.file.menu entryconfig  0 -state active
    .main_menu.file.menu entryconfig  1 -state active
    .main_menu.file.menu entryconfig  2 -state active
    .main_menu.file.menu entryconfig  3 -state active
    .main_menu.file.menu entryconfig  5 -state active
    .main_menu.file.menu entryconfig  6 -state active
    .main_menu.file.menu entryconfig  7 -state active
#    .main_menu.file.menu entryconfig  9 -state active
#    .main_menu.file.menu entryconfig 11 -state active

    .main_menu.help.menu entryconfig  1 -state active
    .main_menu.help.menu entryconfig  2 -state active
    .main_menu.help.menu entryconfig  3 -state active

    ###################################################################
    # Bind all the actions after SCIRun has loaded everything...
    redrawMinicanvas
    bind $minicanvas <B1-Motion> "updateCanvases %x %y"
    bind $minicanvas <1> "updateCanvases %x %y"
    bind $minicanvas <Configure> "redrawMinicanvas"
    $maincanvas bind bgRect <3> "modulesMenu 0 %x %y"
    $maincanvas bind bgRect <1> "startBox $maincanvas %X %Y 0"
    $maincanvas bind bgRect <Control-Button-1> "startBox $maincanvas %X %Y 1"
    $maincanvas bind bgRect <B1-Motion> "makeBox $maincanvas %X %Y"
    $maincanvas bind bgRect <ButtonRelease-1> "$maincanvas delete tempbox"

    # Canvas up-down bound to mouse scroll wheel
    bind . <ButtonPress-5>  "canvasScroll $maincanvas 0.0 0.01"
    bind . <ButtonPress-4>  "canvasScroll $maincanvas 0.0 -0.01"
    
    # Pass in an optional multiplier, since I don't know how to manipulate %D here
    # The %D is typically 120 or -120 on Windows
    bind . <MouseWheel>     "canvasScroll $maincanvas 0.0 %D -.0001"
    bind . <Shift-MouseWheel>     "canvasScroll $maincanvas %D 0.0 -.0001"
    
    # Canvas movement on arrow keys press
    bind . <KeyPress-Down>  "canvasScroll $maincanvas 0.0 0.01"
    bind . <KeyPress-Up>    "canvasScroll $maincanvas 0.0 -0.01"
    bind . <KeyPress-Left>  "canvasScroll $maincanvas -0.01 0.0"
    bind . <KeyPress-Right> "canvasScroll $maincanvas 0.01 0.0" 
    bind . <Destroy> {if {"%W"=="."} {exit 1}} 
    bind . <Control-d> "moduleDestroySelected"
    bind . <Control-l> "ClearCanvas"
    bind . <Control-z> "undo"
    bind . <Control-a> "selectAll"
    bind . <Control-e> "netedit scheduleall"
    bind . <Control-y> "redo"
    bind . <Control-o> "popupLoadMenu"
    bind . <Control-s> "popupSaveMenu"
    bind all <Control-q> "NiceQuit"
}

proc modulesMenu { subnet x y } {
    global mouseX mouseY Subnet
    set mouseX $x
    set mouseY $y
    set canvas $Subnet(Subnet${subnet}_canvas)
    createModulesMenu $canvas.modulesMenu $subnet
    tk_popup $canvas.modulesMenu [expr $x + [winfo rootx $canvas]] \
	[expr $y + [winfo rooty $canvas]]
}

proc redrawMinicanvas {} {
    global SCALEX SCALEY minicanvas maincanvas mainCanvasWidth mainCanvasHeight
    set w [expr [winfo width $minicanvas]-2]
    set h [expr [winfo height $minicanvas]-2]
    set w [expr ($w<=0)?1:$w]
    set h [expr ($h<=0)?1:$h]
    set SCALEX [expr $mainCanvasWidth/$w]
    set SCALEY [expr $mainCanvasHeight/$h]
    updateViewAreaBox
    global Subnet
    set connections ""
    $minicanvas raise module
    foreach module $Subnet(Subnet0_Modules) {
      set coords [scalePath [$maincanvas bbox $module]]
      after 1 $minicanvas coords $module $coords	
      eval lappend connections $Subnet(${module}_connections)
    }
    after 1 drawConnections \{[lsort -unique $connections]\}
    $minicanvas raise "viewAreaBox"
    
    #puts "w:$w h:$h"
    set SCIRUN_SRCDIR [netedit getenv SCIRUN_SRCDIR]
    set lockedimg "$SCIRUN_SRCDIR/pixmaps/locked-24.gif"
    set unlockedimg "$SCIRUN_SRCDIR/pixmaps/unlocked-24.gif"
    #create locked state image

#    global network_executing
#    global cimg
#    if { [string length [info commands ::img::lstate]] } {
#      image delete ::img::lstate
#      $minicanvas delete $cimg
#    } 

#    if { $network_executing } {
#      image create photo ::img::lstate -file $lockedimg
#      set cimg [$minicanvas create image $w $h \
#		      -image ::img::lstate -anchor se]
#    } else {
#      image create photo ::img::lstate -file $unlockedimg
#      set cimg [$minicanvas create image $w $h \
#		      -image ::img::lstate  -anchor se]
#    }

    raise .frame.editor.canvas.mini
}

proc updateViewAreaBox {} {
    global minicanvas maincanvas
    set w [expr [winfo width $minicanvas]-2]
    set h [expr [winfo height $minicanvas]-2]
    set ulx [expr [lindex [$maincanvas xview] 0] * $w]
    set lrx [expr [lindex [$maincanvas xview] 1] * $w]
    set uly [expr [lindex [$maincanvas yview] 0] * $h]
    set lry [expr [lindex [$maincanvas yview] 1] * $h]
    $minicanvas coords viewAreaBox $ulx $uly $lrx $lry

    raise .frame.editor.canvas.mini   
}

proc updateCanvases { x y } {
    global miniCanvasWidth miniCanvasHeight maincanvas minicanvas

    # if the user clicks on the lock / unlock icon, just ignore it
    # if {[expr [winfo width $minicanvas] - $x] < 32 && 
    #    [expr [winfo height $minicanvas] - $y] < 32} return

    set x [expr $x/([winfo width $minicanvas]-2.0)]
    set y [expr $y/([winfo height $minicanvas]-2.0)]
    set xview [$maincanvas xview]
    set yview [$maincanvas yview]
    set x [expr $x-([lindex $xview 1]-[lindex $xview 0])/2]
    set y [expr $y-([lindex $yview 1]-[lindex $yview 0])/2]
    $maincanvas xview moveto $x
    $maincanvas yview moveto $y
    updateViewAreaBox
}

proc createPackageMenu {index} {

    global Color ModuleMenu ModuleMenuHighLight ModuleIPorts ModuleOPorts    
    global ModuleVIPorts ModuleVOPorts time_font
       
    set package [lindex [netedit packageNames] $index]
    set packageToken [join "menu_${package}" ""]
    set ModuleMenu($packageToken) $package
    lappend ModuleMenu(packages) $packageToken
    
    
    foreach category [netedit categoryVNames $package] {	
      set categoryToken [join "${packageToken}_${category}" ""]
      set ModuleMenu($categoryToken) $category
      lappend ModuleMenu(${packageToken}_categories) $categoryToken
      foreach module [netedit moduleVNames $package $category] {
        set moduleToken [join "${categoryToken}_${module}" ""]
        set ModuleMenu($moduleToken) $module
        set ModuleMenuHighLight($moduleToken) 0
        lappend ModuleMenu(${packageToken}_${categoryToken}_modules) $moduleToken
        set "ModuleVIPorts(${package} ${category} ${module})" \
        [netedit module_iport_datatypes $package $category $module]
        set "ModuleVOPorts(${package} ${category} ${module})" \
        [netedit module_oport_datatypes $package $category $module]
      }
    }

    foreach category [netedit categoryNames $package] {	
      foreach module [netedit moduleNames $package $category] {
        set "ModuleIPorts(${package} ${category} ${module})" \
        [netedit module_iport_datatypes $package $category $module]
        set "ModuleOPorts(${package} ${category} ${module})" \
        [netedit module_oport_datatypes $package $category $module]
      }
    }

    set pack $packageToken
    # Add the cascade button and menu for the package to the menu bar

    .main_menu.module.menu configure -disabledforeground $Color(MenuHeadingForeGround) 
    if { $index == 1} {
      .main_menu.module.menu add command -label "$ModuleMenu($pack)" -state disabled  
    } else {
      .main_menu.module.menu add command -label "$ModuleMenu($pack)" -state disabled
    }
    
    foreach cat $ModuleMenu(${pack}_categories) {
      # Add the category to the menu bar menu

      .main_menu.module.menu add cascade -label "$ModuleMenu($cat)" \
        -menu .main_menu.module.menu.$cat 

      menu .main_menu.module.menu.$cat -tearoff true\
        -bd 2 -activeborderwidth 0 \
        -activebackground $Color(MenuSelectBackGround)  \
        -activeforeground $Color(MenuSelectForeGround) \
        -background $Color(MenuBackGround) \
        -foreground $Color(MenuForeGround)

      foreach mod $ModuleMenu(${pack}_${cat}_modules) {
          .main_menu.module.menu.$cat add command \
            -label "$ModuleMenu($mod)" \
            -command "addModule \"$ModuleMenu($pack)\" \"$ModuleMenu($cat)\" \"$ModuleMenu($mod)\""
      }
    } 

    update idletasks
}

proc UpdateHighLightsPackageMenu {} {
    global Color ModuleMenu ModuleMenuHighLight
    
    HighLightModules
    
    set packages [netedit packageNames]
    set cindex 0
    
    foreach pack $packages {
    
    set pack [join "menu_${pack}" ""]    
      incr cindex
    
      foreach cat $ModuleMenu(${pack}_categories) {
        # Add the category to the menu bar menu

        set mindex 1
        set foundhl 0

        foreach mod $ModuleMenu(${pack}_${cat}_modules) {
        
          if { $ModuleMenuHighLight($mod) } {
            .main_menu.module.menu.$cat entryconfigure $mindex  -background $Color(MenuBackGroundHL) -foreground $Color(MenuForeGroundHL) -activebackground $Color(MenuSelectBackGroundHL) -activeforeground $Color(MenuSelectForeGroundHL)
            set foundhl 1
          } else {
            .main_menu.module.menu.$cat entryconfigure $mindex -background $Color(MenuBackGround) -foreground $Color(MenuForeGround) -activebackground $Color(MenuSelectBackGround) -activeforeground $Color(MenuSelectForeGround)
          }
          incr mindex
        }
        if {[expr $foundhl == 1]} {
          .main_menu.module.menu entryconfigure $cindex -background $Color(MenuBackGroundHL) -foreground $Color(MenuForeGroundHL) -activebackground $Color(MenuSelectBackGroundHL) -activeforeground $Color(MenuSelectForeGroundHL)
        } else {
          .main_menu.module.menu entryconfigure $cindex -background $Color(MenuBackGround) -foreground $Color(MenuForeGround) -activebackground $Color(MenuSelectBackGround) -activeforeground $Color(MenuSelectForeGround)
        }
        incr cindex
      } 
    }
}

proc HighLightModules {} {
  global Color ModuleMenu ModuleMenuHighLight ModuleFilter
  
  set names [array names ModuleMenuHighLight]
  if {![string length $ModuleFilter]} {
    foreach name $names {
      set ModuleMenuHighLight($name) 0
    }
   return;
  }
  
  foreach name $names {
    set ModuleMenuHighLight($name) 1
  }
  foreach patt $ModuleFilter {
    foreach name $names {
      set modname [set ModuleMenu($name)]
      if {![string match -nocase "*$patt*" $modname]} {
        set ModuleMenuHighLight($name) 0
      }
    }
  }
}


# createModulesMenu is called when the user right-clicks on the
# canvas.  It presents them with a menu of all modules.  Selecting
# a module name will create it at the clicked location
proc createModulesMenu { menu subnet } {
    global ModuleMenu Subnet Color
    # return if there is no information to put in menu
    if ![info exists ModuleMenu] return
    # destroy the old menu
    #    if [winfo exists $menu] {	
    #	destroy $menu
    #    }
    # create a new menu
    if { ![winfo exists $menu] } {	
      menu $menu -bd 2 -activeborderwidth 0 \
        -activebackground $Color(MenuSelectBackGround)  \
        -activeforeground $Color(MenuSelectForeGround) \
        -background $Color(MenuBackGround) \
        -foreground $Color(MenuForeGround)
    }
    $menu delete 0 end
    $menu configure -tearoff false -disabledforeground black

    foreach pack $ModuleMenu(packages) {
      # Add a menu separator if this package isn't the first one
      if { [$menu index end] != "none" } {
          $menu add separator 
      }
      # Add a label for the Package name
      $menu add command -label "$ModuleMenu($pack)" -state disabled
      foreach cat $ModuleMenu(${pack}_categories) {
        # Add the category to the right-button menu
        $menu add cascade -label "  $ModuleMenu($cat)" -menu $menu.$cat
        if { ![winfo exists $menu.$cat] } {	
          menu $menu.$cat -tearoff false -bd 2 -activeborderwidth 0 \
            -activebackground $Color(MenuSelectBackGround)  \
            -activeforeground $Color(MenuSelectForeGround) \
            -background $Color(MenuBackGround) \
            -foreground $Color(MenuForeGround)
        }
        $menu.$cat delete 0 end

        foreach mod $ModuleMenu(${pack}_${cat}_modules) {
          $menu.$cat add command -label "$ModuleMenu($mod)" \
            -command "addModuleAtMouse \"$ModuleMenu($pack)\" \"$ModuleMenu($cat)\" \"$ModuleMenu($mod)\" \"$subnet\""
        }
      }
    }
    
   $menu add separator
    $menu add cascade -label "Subnets" -menu $menu.subnet
    if { ![winfo exists $menu.subnet] } {	
      menu $menu.subnet -tearoff false -bd 2 -activeborderwidth 0 \
        -activebackground $Color(MenuSelectBackGround)  \
        -activeforeground $Color(MenuSelectForeGround) \
        -background $Color(MenuBackGround) \
        -foreground $Color(MenuForeGround)
    }

    createSubnetMenu $menu $subnet
	
    $menu configure -postcommand "UpdateHighLightsModulesMenu $menu $subnet"
    update idletasks
}

proc UpdateHighLightsModulesMenu { menu subnet } {
    global Color ModuleMenu ModuleMenuHighLight
    
    HighLightModules
    
    set cindex -1
    
    foreach pack $ModuleMenu(packages) {
      incr cindex
      incr cindex
    
      foreach cat $ModuleMenu(${pack}_categories) {
        # Add the category to the menu bar menu

        set mindex 0
        set foundhl 0

        foreach mod $ModuleMenu(${pack}_${cat}_modules) {
          
          if { $ModuleMenuHighLight($mod) } {
            $menu.$cat entryconfigure $mindex  -background $Color(MenuBackGroundHL) -foreground $Color(MenuForeGroundHL) -activebackground $Color(MenuSelectBackGroundHL) -activeforeground $Color(MenuSelectForeGroundHL)
            set foundhl 1
          } else {
            $menu.$cat entryconfigure $mindex -background $Color(MenuBackGround) -foreground $Color(MenuForeGround) -activebackground $Color(MenuSelectBackGround) -activeforeground $Color(MenuSelectForeGround)
          }
          incr mindex
        }
        
        if {[expr $foundhl == 1]} {
          $menu entryconfigure $cindex -background $Color(MenuBackGroundHL) -foreground $Color(MenuForeGroundHL) -activebackground $Color(MenuSelectBackGroundHL) -activeforeground $Color(MenuSelectForeGroundHL)
        } else {
          $menu entryconfigure $cindex -background $Color(MenuBackGround) -foreground $Color(MenuForeGround) -activebackground $Color(MenuSelectBackGround) -activeforeground $Color(MenuSelectForeGround)
        }
        incr cindex
      } 
    }
}


proc createSubnetMenu { { menu "" } { subnet 0 } } {
    global SubnetScripts time_font
    loadSubnetScriptsFromDisk
    #generateSubnetScriptsFromNetwork

    if { [winfo exists $menu ] } {
      $menu.subnet delete 0 end
    }
    .main_menu.subnet.menu delete 0 end
    set names [lsort -dictionary [array names SubnetScripts *]]

    if { ![llength $names] } {
#      if { [winfo exists $menu ] } {
#          $menu entryconfigure Subnets -state disabled
#      }
#      .main_menu configure Subnets configure -state disabled
    } else {
      if { [winfo exists $menu ] } {
          $menu entryconfigure Subnets -state normal
      }
      .main_menu.subnet configure -state normal
      foreach name $names {
          if { [winfo exists $menu ] } {
            $menu.subnet add command -label "$name" \
            -command "instanceSubnet \"$name\" 0 0 $subnet"
          }
          .main_menu.subnet.menu add command -label "$name" \
            -command "instanceSubnet \"$name\" 0 0 $subnet" 
      }
    }
}
    
proc createToolkitsMenu {} {
    global Color
    .main_menu.toolkits.menu delete 0 end

    set srcdir [netedit getenv "SCIRUN_EXAMPLE_NETS_DIR"]
    set toolkit $srcdir/FwdInvToolbox

    set toolkit_subdirs [glob -directory $toolkit -type d *]
    .main_menu.toolkits.menu add cascade -label "ECG Forward/Inverse Toolkit" -menu .main_menu.toolkits.menu.fwdinvtoolkit
    if { ![winfo exists .main_menu.toolkits.menu.fwdinvtoolkit] } {
      menu .main_menu.toolkits.menu.fwdinvtoolkit -tearoff false -bd 2 \
        -activeborderwidth 0 \
        -activebackground $Color(MenuSelectBackGround) \
        -activeforeground $Color(MenuSelectForeGround) \
        -background $Color(MenuBackGround) \
        -foreground $Color(MenuForeGround)
    }
    set fwdinvtoolkit_menu .main_menu.toolkits.menu.fwdinvtoolkit
    foreach subdir $toolkit_subdirs {
      set fragments [split $subdir "/\\"]
      set subdir_menu [lindex $fragments end]

      if { ![winfo exists $fwdinvtoolkit_menu.$subdir_menu] } {
        $fwdinvtoolkit_menu add cascade -label "$subdir_menu" -menu $fwdinvtoolkit_menu.$subdir_menu

        set toolkit_nets [glob -directory $subdir *.srn]
        menu $fwdinvtoolkit_menu.$subdir_menu -tearoff false -bd 2 \
          -activeborderwidth 0 \
          -activebackground $Color(MenuSelectBackGround) \
          -activeforeground $Color(MenuSelectForeGround) \
          -background $Color(MenuBackGround) \
          -foreground $Color(MenuForeGround)

        foreach net $toolkit_nets {
	  set net_fragments [split $net "/\\."]
          set net_name [lindex $net_fragments end-1]
          $fwdinvtoolkit_menu.$subdir_menu add command -label "$net_name" -command "loadToolkitNet \"$net\""
        }
      }
    }
}

proc networkHasChanged {args} {
    upvar \#0 NetworkChanged changed
#    puts "$changed networkHasChanged [info level [expr [info level]-1]]"
    set changed 1
}

proc addModule { package category module } {
    return [addModuleAtPosition "$package" "$category" "$module" 10 10]
}

proc addModuleAtMouse { pack cat mod subnet_id } {
    global mouseX mouseY Subnet
    set Subnet(Loading) $subnet_id
    set ret [addModuleAtPosition $pack $cat $mod $mouseX $mouseY]
    set Subnet(Loading) 0
    set mouseX 10
    set mouseY 10
    return $ret
}


proc findMovedModulePath { packvar catvar modvar } {
    # Deprecated module translation table.
}	        


proc addModuleAtPosition {package category module { xpos 10 } { ypos 10 } { absolute 0 } { modid "" } } {
    # Look up the real category for a module.  This allows networks to
    # be read in if the modules change categories.
    findMovedModulePath package category module
    set category [netedit getCategoryName $package $category $module]

    # fix for bug #2052, allows addmodule to call undefined procs without error
    set unknown_body [info body unknown]
    proc unknown { args } {}
    
    # default argument is empty, but if C already created the module, 
    # it will pass the id in.
    if { $modid == "" } {
      # Tell the C++ network to create the requested module
      set modid [netedit addmodule "$package" "$category" "$module"]
    }
    # Reset the unknown proc to default behavior
    proc unknown { args } $unknown_body

    # netedit addmodule returns an empty string if the module wasnt created
    if { ![string length $modid] } {
      tk_messageBox -type ok -parent . -icon warning -message \
        "Cannot find the ${package}::${category}::${module} module."
      return
    }    

    networkHasChanged
    global inserting Subnet
    set canvas $Subnet(Subnet$Subnet(Loading)_canvas)
    set Subnet($modid) $Subnet(Loading)
    set Subnet(${modid}_connections) ""
    lappend Subnet(Subnet$Subnet(Loading)_Modules) $modid

    set className [join "${package}_${category}_${module}" ""]
    # Create the itcl object
    if {[catch "$className $modid" exception]} {
      # Use generic module
      if {$exception != "invalid command name \"$className\""} {
          puts "Error instantiating iTcl class for module:\n$exception";
      }
      Module $modid -name "$module" -has_ui 0
    }

    # compute position if we're inserting the net to the right    
    if { $inserting } {
      global insertOffset
      set xpos [expr $xpos+[lindex $insertOffset 0]]
      set ypos [expr $ypos+[lindex $insertOffset 1]]
    } else { ;# create the module relative to current screen position
      set xpos [expr $xpos+[$canvas canvasx 0]]
      set ypos [expr $ypos+[$canvas canvasy 0]]
    }
    set absolute [expr $absolute || $inserting]
    $modid make_icon $xpos $ypos $absolute
    update idletasks

    raise .frame.editor.canvas.mini
    return $modid    
}


proc addModuleAtAbsolutePosition {package category module { xpos 10 } { ypos 10 } { modid "" } } {
    # Look up the real category for a module.  This allows networks to
    # be read in if the modules change categories.
    findMovedModulePath package category module
    set category [netedit getCategoryName $package $category $module]

    # fix for bug #2052, allows addmodule to call undefined procs without error
    set unknown_body [info body unknown]
    proc unknown { args } {}
    
    # default argument is empty, but if C already created the module, 
    # it will pass the id in.
    if { $modid == "" } {
      # Tell the C++ network to create the requested module
      set modid [netedit addmodule "$package" "$category" "$module"]
    }
    # Reset the unknown proc to default behavior
    proc unknown { args } $unknown_body

    # netedit addmodule returns an empty string if the module wasnt created
    if { ![string length $modid] } {
      tk_messageBox -type ok -parent . -icon warning -message \
          "Cannot find the ${package}::${category}::${module} module."
      return
    }    

    networkHasChanged
    global inserting Subnet
    set canvas $Subnet(Subnet$Subnet(Loading)_canvas)
    set Subnet($modid) $Subnet(Loading)
    set Subnet(${modid}_connections) ""
    lappend Subnet(Subnet$Subnet(Loading)_Modules) $modid

    set className [join "${package}_${category}_${module}" ""]
    # Create the itcl object
    if {[catch "$className $modid" exception]} {
      # Use generic module
      if {$exception != "invalid command name \"$className\""} {
          puts "Error instantiating iTcl class for module:\n$exception";
      }
      Module $modid -name "$module" -has_ui 0
    }

    # compute position if we're inserting the net to the right    
    if { $inserting } {
      global insertOffset
      set xpos [expr $xpos+[lindex $insertOffset 0]]
      set ypos [expr $ypos+[lindex $insertOffset 1]]
    } else { ;# create the module relative to current screen position
    }
    $modid make_icon $xpos $ypos 1
    update idletasks

    raise .frame.editor.canvas.mini    
    return $modid
}


proc append_srn_filename {name} {

    set ext_ind [expr [string length $name] - 4]
    set ext [string range $name $ext_ind end]
    
    if { $ext == ".net" } {
      set name [string range $name 0 $ext_ind]srn
      createSciDialog -warning -title "Save Warning" -button1 "Ok"\
        -message "SCIRun no longer saves .net files.\nSaving $name instead."
      set ext ".srn"
    } 
    
    if { $ext != ".srn" } {
      set name $name.srn
    } 
    return $name
}


proc popupSaveMenu {} {
    global netedit_savefile NetworkChanged
    if { $netedit_savefile != "" } {
      # We know the name of the savefile, dont ask for name, just save it
      # make sure we only save .srn files
      set netedit_savefile [append_srn_filename $netedit_savefile]
      wm title . "SCIRun ([lindex [file split $netedit_savefile] end])"
      writeNetwork $netedit_savefile
      netedit setenv "SCIRUN_NETFILE" $netedit_savefile
      set NetworkChanged 0
    } else { ;# Otherwise, ask the user for the name to save as
      popupSaveAsMenu
    }
}

proc popupSaveAsMenu {} {
    set types {
      {{SCIRun Net} {.srn} }
    } 

    global netedit_savefile NetworkChanged

    # determine initialdir based on current $netedit_savefile
    set dirs [file split "$netedit_savefile"]
    set initialdir [pwd]
    set initialfile ""

    if {[llength $dirs] > 1} {
      set initialdir ""
      set size [expr [llength $dirs] - 1]
      for {set i 0} {$i<$size} {incr i} {
          set initialdir [file join $initialdir [lindex $dirs $i]]
      }
      set initialfile [lindex $dirs $size]
    }

    global save_file netedit_savefile filetype_savefile
    
    set save_file 0;
    set netedit_savefile $initialfile
    set filetype_savefile ".srn"

    set w .save_network_dialog
    if {![winfo exists $w]} {
      sci_toplevel $w -class TkFDialog
      
      makeSaveFilebox \
          -parent $w \
          -filevar netedit_savefile \
          -cancel "global save_file; set save_file 0; wm withdraw $w" \
          -commandname "Save" \
          -command "global save_file; set save_file 1; wm withdraw $w" \
          -title "Save Network file" \
          -filetypes $types \
          -initialfile $initialfile \
          -initialdir $initialdir \
          -formatvar filetype_savefile \
          -defaultextension ".srn" \
          -formats {None} \
          -donotresize 1
    } else {
      wm deiconify $w
    }

    tkwait variable save_file
    
    if {$save_file == 1} {
      if { $netedit_savefile != "" } {
        # make sure we only save .srn files
        set netedit_savefile [append_srn_filename $netedit_savefile]
        wm title . "SCIRun ([lindex [file split $netedit_savefile] end])"
        writeNetwork $netedit_savefile
        set NetworkChanged 0
        netedit setenv "SCIRUN_NETFILE" $netedit_savefile
      }
    }
}

proc popupInsertMenu { {subnet 0} } {
    global inserting insertOffset Subnet NetworkChanged
    global mainCanvasWidth mainCanvasHeight
    
    #get the net to be inserted
    set types {
      {{SCIRun Net} {.srn} }
      {{old SCIRun Net} {.net} }
    } 

    global insert_file netedit_insertfile filetype_insertfile
    set w .insert_network_dialog

    set initialdir [netedit getenv SCIRUN_NET_DIR] 
    set insert_file 0
    set netedit_insertfile ""
    set filetype_insertfile ".srn"   

    if {![winfo exists $w]} {
      sci_toplevel $w -class TkFDialog
      makeOpenFilebox \
        -parent $w \
        -filevar netedit_insertfile \
        -cancel "global insert_file; set insert_file 0; wm withdraw $w" \
        -commandname "Insert" \
        -command "global insert_file; set insert_file 1; wm withdraw $w" \
        -title "Insert Network file" \
        -filetypes $types \
        -initialdir $initialdir \
        -defaultextension ".srn" \
        -donotresize 1        
    } else {
      wm deiconify $w
    }

    tkwait variable insert_file
    
    if {$insert_file != 1} return
 
#    set netedit_insertfile [tk_getOpenFile -filetypes $types -initialdir [netedit getenv SCIRUN_NET_DIR] ]
    if { [check_filename $netedit_insertfile] == "invalid" } {
      set netedit_insertfile ""
      return
    }
    if { ![file exists $netedit_insertfile]} { 
      return
    }
    
    set canvas $Subnet(Subnet${subnet}_canvas)    
    # get the bbox for the net being inserted by
    # parsing netedit_openfile for bbox 
    set fchannel [open $netedit_insertfile]
    set curr_line ""
    set curr_line [gets $fchannel]
    while { ![eof $fchannel] } {
    if { [string match "set bbox*" $curr_line] } {
	    eval $curr_line
	    break
    }
    set curr_line [gets $fchannel]
    set val [string match "set bbox*" $curr_line]

    }
    set viewBox "0 0 [winfo width $canvas] [winfo width $canvas]"
    if { ![info exists bbox] || [llength $bbox] != 4 } {
      set bbox $viewBox
    }
    set w [expr [lindex $bbox 2] - [lindex $bbox 0]]
    set h [expr [lindex $bbox 3] - [lindex $bbox 1]]
    set oldbox $bbox
    set moveBoxX "[expr $w/2] 0 [expr $w/2] 0"
    set moveBoxY "0 [expr $h/2] 0 [expr $h/2]"
    set done 0
    set bbox [list 0 0 $w $h] ;# start inserting in upper left corner
    while {!$done} {
      set done 1
      set modules [eval $canvas find overlapping $bbox]
      foreach modid $modules {
        if { [lsearch [$canvas gettags $modid] module] != -1 } {
          set overlap [clipBBoxes [compute_bbox $canvas $modid] $bbox]
          if ![string equal $overlap "0 0 0 0"] {
            set done 0
            break
          }
        }
      }
      if {!$done} {
        # move the insert position left by half a screen
        for {set i 0} {$i < 4} {incr i} {
          set bbox [lreplace $bbox $i $i \
			      [expr [lindex $moveBoxX $i]+[lindex $bbox $i]]]
        }
        if {[lindex $bbox 2] > $mainCanvasWidth} {
          set bbox [lreplace $bbox 2 2 \
          [expr [lindex $bbox 2] -[lindex $bbox 0]]]
          set bbox [lreplace $bbox 0 0 0]
          for {set i 0} {$i < 4} {incr i} {
            set bbox [lreplace $bbox $i $i \
			       [expr [lindex $moveBoxY $i]+[lindex $bbox $i]]]
          }
          if {[lindex $bbox 3] > $mainCanvasHeight} {
            set bbox [list 50 50 [expr 50+$w] [expr 50+$h]]
            set done 1
          }
        }
      }
    }

    set insertOffset [list [expr [lindex $bbox 0]-[lindex $oldbox 0]] \
			  [expr [lindex $bbox 1]-[lindex $oldbox 1]]]
    $canvas xview moveto [expr [lindex $bbox 0]/$mainCanvasWidth-0.01]
    $canvas yview moveto [expr [lindex $bbox 1]/$mainCanvasHeight-0.01]
    set preLoadModules $Subnet(Subnet${subnet}_Modules)
    set inserting 1
    if {[string match *.net $netedit_insertfile]} {
      loadnet $netedit_insertfile
    } else {
      global netedit_savefile
      set tmp $netedit_savefile
      uplevel \#0 netedit load_srn $netedit_insertfile
      set netedit_savefile $tmp
    }
    set inserting 0
    unselectAll
    foreach module $Subnet(Subnet${subnet}_Modules) {
      if { [lsearch $preLoadModules $module] == -1 } {
        $module addSelected
      }
    }
    
    set NetworkChanged 1
}

proc subnet_bbox { subnet { cheat 1} } {
    global Subnet
    return [compute_bbox $Subnet(Subnet${subnet}_canvas) \
		$Subnet(Subnet${subnet}_Modules) 1]
}

proc compute_bbox { canvas { items "" } { cheat 0 } } {
    set canvasbounds [$canvas cget -scrollregion]
    set maxx [lindex $canvasbounds 0]
    set maxy [lindex $canvasbounds 1]
    set minx [lindex $canvasbounds 2]
    set miny [lindex $canvasbounds 3]
    global CurrentlySelectedModules
    if { $items == ""} { set items $CurrentlySelectedModules }
    if { ![llength $items] } { return [list 0 0 0 0] }
    foreach item $items {
	set bbox [$canvas bbox $item]
	if { $cheat && [lsearch [$canvas gettags $item] module] != -1 } {
	    if { [expr [lindex $bbox 2] - [lindex $bbox 0] < 10]  } {
		set bbox [lreplace $bbox 2 2 [expr [lindex $bbox 0]+180]]
	    }
	    if { [expr [lindex $bbox 3] - [lindex $bbox 1] < 10]  } {
		set bbox [lreplace $bbox 3 3 [expr [lindex $bbox 1]+80]]
	    }
	}
	if { [lindex $bbox 0] <= $minx} { set minx [lindex $bbox 0] }
	if { [lindex $bbox 1] <= $miny} { set miny [lindex $bbox 1] }
	if { [lindex $bbox 2] >  $maxx} { set maxx [lindex $bbox 2] }
	if { [lindex $bbox 3] >  $maxy} { set maxy [lindex $bbox 3] }
    }
    return [list $minx $miny $maxx $maxy]
}


proc check_filename {name} {

	if {$name == ""} {
		return "invalid"
	}

    set ext_ind [expr [string length $name] - 4]
    set ext [string range $name $ext_ind end]
    
    if { $ext != ".net" && $ext != ".srn"} {
		set name [string range $name 0 $ext_ind]srn
		set msg "Valid net files end with .srn (or .net prior to v1.25.2)"
		createSciDialog -warning -title "Save Warning" -button1 "Ok"\
			-message $msg
		return "invalid"
    } 
    return "valid"
}

proc update_network_editor_title {filename} {
    # remove full path stuff
    set filename [lindex [file split $filename] end]
    wm title . "SCIRun ($filename)"
}

proc popupLoadMenu {} {
    global NetworkChanged
    if $NetworkChanged {
      set result [createSciDialog -warning -title "Save Network?" \
		    -button1 "No" -button2 "Cancel" -button3 "Yes" \
		    -message "Your network has not been saved.\n\nWould you like to save before loading a new one?"]    
    
      if {$result == 3} { popupSaveMenu }
      if {$result == 2} { return }
    }

    set types {
      {{SCIRun Net} {.srn} }
      {{old SCIRun Net} {.net} }
    } 

    global open_file netedit_openfile
    set initialdir [netedit getenv SCIRUN_NET_DIR] 
    set open_file 0
    set netedit_openfile ""

    set w .open_network_dialog
    if {![winfo exists $w]} {
      sci_toplevel $w -class TkFDialog
      makeOpenFilebox \
          -parent $w \
          -filevar netedit_openfile \
          -cancel "global open_file; set open_file 0; wm withdraw $w" \
          -commandname "Open" \
          -command "global open_file; set open_file 1; wm withdraw $w" \
          -title "Open Network file" \
          -filetypes $types \
          -initialdir $initialdir \
          -defaultextension ".srn" \
          -donotresize 1
    } else {
      wm deiconify $w
    }

    tkwait variable open_file
    
    if {$open_file != 1} return
 
    puts " $netedit_openfile"
 
    if { [check_filename $netedit_openfile] == "invalid" } {
      set netedit_openfile ""
      return
    }

    if { ![file exists $netedit_openfile]} { 
      return
    }

    #dont ask user before clearing canvas
    ClearCanvas 0
    set inserting 0
    if {[string match *.srn $netedit_openfile]} {
      # compensate for spaces in the filename (windows)
      after 500 uplevel \#0 netedit load_srn \{\{$netedit_openfile\}\}
    } else {
      loadnet $netedit_openfile 
    }
    wm title . "SCIRun ([lindex [file split $netedit_openfile] end])"
}

proc loadToolkitNet { net } {
  # basically copied from popupLoadMenu
  global NetworkChanged
  if $NetworkChanged {
    set result [createSciDialog -warning -title "Save Network?" \
		  -button1 "No" -button2 "Cancel" -button3 "Yes" \
		  -message "Your network has not been saved.\n\nWould you like to save before loading a new one?"]    
    
    if {$result == 3} { popupSaveMenu }
    if {$result == 2} { return }
  }
  puts " $net"
 
  if { [check_filename $net] == "invalid" } {
    return
  }

  if { ![file exists $net]} { 
    return
  }

  # don't ask user before clearing canvas
  ClearCanvas 0
  set inserting 0
  if {[string match *.srn $net]} {
    # compensate for spaces in the filename (windows)
    after 500 uplevel \#0 netedit load_srn \{\{$net\}\}
  } else {
    loadnet $net
  }
  wm title . "SCIRun ([lindex [file split $net] end])"
}

proc ClearCanvas { { confirm 1 } { subnet 0 } } {
    # destroy all modules
    global NetworkChanged
    if { !$NetworkChanged } { set confirm 0 }
    set do_clear 3

    if { $confirm } {
      set message [list "Your network has not been saved." \
			 "All Modules and connections will be deleted." \
			 "Really clear?"]
      set do_clear [createSciDialog -warning -title "Warning" \
		    -button1 "No" -button3 "Yes" \
		    -message [join $message "\n\n"]]
    }

  if { $do_clear == 3 } {
    global Subnet
    foreach module $Subnet(Subnet${subnet}_Modules) {
	    if { [string first Render_Viewer $module] != -1 } {
        moduleDestroy $module
        }
      }
      foreach module $Subnet(Subnet${subnet}_Modules) {
        moduleDestroy $module
      }

      wm title . "SCIRun" ;# Reset Main Window Title
      setGlobal netedit_savefile ""
      setGlobal CurrentlySelectedModules ""
      setGlobal NetworkChanged 0
	  netedit setenv "SCIRUN_NETFILE" ""
    
      set Subnet(num) 0
    }
    
    raise .frame.editor.canvas.mini
}

proc UpdateCanvas { {subnet 0 } } {

    global Subnet
    foreach module $Subnet(Subnet${subnet}_Modules) {
      if [isaSubnetIcon $module] {
        UpdateCanvas $Subnet(${module}_num)
      }
      $module update_icon
    }

    set connections ""
    foreach module $Subnet(Subnet${subnet}_Modules) {
      eval lappend connections $Subnet(${module}_connections)
    }
    
    if {[expr [string compare $connections ""] == 0]} {
        after 1 drawConnections \{[lsort -unique $connections]\}
    }
    
    if {[expr $subnet == 0]} {
      raise .frame.editor.canvas.mini
    }
}

proc NiceQuit {} {
    global NetworkChanged netedit_savefile
    if { $NetworkChanged && ![envBool SCIRUN_FAST_QUIT] } {
      set result [createSciDialog -warning -title "Quit?" \
		    -button1 "Don't Save" -button2 "Cancel" -button3 "Save" \
		    -message "Your session has not been saved.\nWould you like to save before exiting?"]
      switch -- $result { 
        "-1" return
        "2" return
        "3" {
        if { [winfo exists .standalone] } {
              app save_session
          } else {
              puts -nonewline "Saving $netedit_savefile..."
              popupSaveMenu
          }
        }
      }
    }
    set geom [open ~/.scirun.geom w]
    puts $geom [wm geom .]
    close $geom
    puts "Goodbye!"
    netedit quit
}


proc UpdateGraphicsDriversDialog {} {
    createSciDialog -error -title "Please Update Graphics Card Drivers" \
      -button1 "OK" \
      -message "The graphics card drivers on this computer are out of date.\nSCIRun requires OpenGL version 2.0 or higher.\nPlease download an updated graphics card driver for this machine.\nAn up to date driver can be generally found at the website of graphics card manufacturer.\n"
}


proc initInfo { {subnet_number 0} } {
    global Subnet
    if { ![info exists Subnet(Subnet${subnet_number}_userName)] } {
	set Subnet(Subnet${subnet_number}_userName) [netedit getenv LOGNAME] 
	if { $Subnet(Subnet${subnet_number}_userName) == "" } {
	    set Subnet(Subnet${subnet_number}_userName) [netedit getenv USER] 
	}
	if { $Subnet(Subnet${subnet_number}_userName) == "" } {
	    set Subnet(Subnet${subnet_number}_userName) Unknown
	}
    }
    if { ![info exists Subnet(Subnet${subnet_number}_creationDate)] } {
	set Subnet(Subnet${subnet_number}_creationDate) \
	    [clock format [clock seconds] -format "%a %b %d %Y"]
    }
    if { ![info exists Subnet(Subnet${subnet_number}_creationTime)] } {
	set Subnet(Subnet${subnet_number}_creationTime) \
	    [clock format [clock seconds] -format "%H:%M:%S"]
    }
    if { ![info exists Subnet(Subnet${subnet_number}_runDate)] } {
	set Subnet(Subnet${subnet_number}_runDate) ""
    }
    if { ![info exists Subnet(Subnet${subnet_number}_runTime)] } {
	set Subnet(Subnet${subnet_number}_runTime) ""
    }

    if { ![info exists Subnet(Subnet${subnet_number}_notes)] } { 
	set Subnet(Subnet${subnet_number}_notes) "" 
    }
		
		if { ![info exists Subnet(Subnet${subnet_number}_relfilenames)]} {
			set Subnet(Subnet${subnet_number}_relfilenames) 0
		}

		if { ![info exists Subnet(Subnet${subnet_number}_netversion)]} {
			set Subnet(Subnet${subnet_number}_netversion) 0
		}
}

proc updateRunDateAndTime { subnet_number } {
    global Subnet
    set Subnet(Subnet${subnet_number}_runDate) \
	[clock format [clock seconds] -format "%a %b %d %Y"]
    set Subnet(Subnet${subnet_number}_runTime) \
	[clock format [clock seconds] -format "%H:%M:%S"]
}    


proc updateCreationDateAndTime { subnet_number } {
    global Subnet
    set Subnet(Subnet${subnet_number}_creationDate) \
	[clock format [clock seconds] -format "%a %b %d %Y"]
    set Subnet(Subnet${subnet_number}_creationTime) \
	[clock format [clock seconds] -format "%H:%M:%S"]
}    

proc getNetVersion {} {
  global Subnet
  return [set Subnet(Subnet0_netversion)]
}

# This proc was added by Mohamed Dekhil to save some info about the net
proc popupInfoMenu { {subnet_num 0 } } {
    global Subnet 
    initInfo $subnet_num
    
    lappend backupDateTimeNotes \
      $Subnet(Subnet${subnet_num}_userName) \
      $Subnet(Subnet${subnet_num}_runDate) \
      $Subnet(Subnet${subnet_num}_runTime) \
      $Subnet(Subnet${subnet_num}_creationDate) \
      $Subnet(Subnet${subnet_num}_creationTime) \
      $Subnet(Subnet${subnet_num}_notes) \
      $Subnet(Subnet${subnet_num}_relfilenames) \
      $Subnet(Subnet${subnet_num}_netversion) \

    set w .netedit_info$subnet_num
        
    if {[winfo exists $w]} {
			raise $w
			return
    }

    toplevel $w
    wm title $w "$Subnet(Subnet${subnet_num}_Name) Properties"

		iwidgets::tabnotebook $w.tab -tabpos n  -width 600 -height 600
		
		pack $w.tab -fill both -expand yes
		set infoframe [$w.tab add -label "Information" ]
		# set saveframe [$w.tab add -label "Save Options" ] 
		
		$w.tab select 0

		# disable for now until we resolve problems with srn files
		# frame $saveframe.fname
		# checkbutton $saveframe.fname.relativefiles -variable Subnet(Subnet${subnet_num}_relfilenames) -text "Save all filenames relative to the network location"
		
		# pack $saveframe.fname -side top -padx 1 -pady 1 -ipadx 2 -ipady 2 -fill x
		# pack $saveframe.fname.relativefiles -side left

    frame $infoframe.fname
    label $infoframe.fname.lname -text "User: " -padx 3 -pady 3
    entry $infoframe.fname.ename -width 50 -relief sunken -bd 2 \
      -textvariable Subnet(Subnet${subnet_num}_userName)
    pack $infoframe.fname.lname $infoframe.fname.ename -side left

    set pre [expr $subnet_num?"Sub-":""]
    frame $infoframe.cdt
    label $infoframe.cdt.label -text "${pre}Network Created:"
    label $infoframe.cdt.ldate -text "   Date: " -padx 3 -pady 3 
    entry $infoframe.cdt.edate -width 20 -relief sunken -bd 2 \
      -textvariable Subnet(Subnet${subnet_num}_creationDate)

    label $infoframe.cdt.ltime -text "   Time: " -padx 5 -pady 3 
    entry $infoframe.cdt.etime -width 10 -relief sunken -bd 2 \
      -textvariable Subnet(Subnet${subnet_num}_creationTime)

    button $infoframe.cdt.reset -text "Reset" \
      -command "updateCreationDateAndTime $subnet_num"
    pack $infoframe.cdt.label  -side left -fill x
    pack $infoframe.cdt.reset $infoframe.cdt.etime $infoframe.cdt.ltime $infoframe.cdt.edate $infoframe.cdt.ldate   -side right -padx 5


    frame $infoframe.fdt
    label $infoframe.fdt.label -text "${pre}Network Executed:"
    label $infoframe.fdt.ldate -text "   Date: " -padx 3 -pady 3 
    entry $infoframe.fdt.edate -width 20 -relief sunken -bd 2 \
      -textvariable Subnet(Subnet${subnet_num}_runDate)

    label $infoframe.fdt.ltime -text "   Time: " -padx 5 -pady 3 
    entry $infoframe.fdt.etime -width 10 -relief sunken -bd 2 \
      -textvariable Subnet(Subnet${subnet_num}_runTime)

    button $infoframe.fdt.reset -text "Reset" \
      -command "updateRunDateAndTime $subnet_num"

    pack $infoframe.fdt.label  -side left -fill x
    pack $infoframe.fdt.reset $infoframe.fdt.etime $infoframe.fdt.ltime $infoframe.fdt.edate $infoframe.fdt.ldate   -side right -padx 5

    frame $infoframe.nv
    label $infoframe.nv.label -text "Network version: "
    entry $infoframe.nv.entry -width 20 -relief sunken -bd 2 \
      -textvariable Subnet(Subnet${subnet_num}_netversion)

    pack $infoframe.nv.entry $infoframe.nv.label -side right -padx 5


    frame $infoframe.fnotes -relief groove
    frame $infoframe.fnotes.top
    frame $infoframe.fnotes.bot
    label $infoframe.fnotes.top.lnotes -text "Notes:" -padx 2 -pady 5 
    text $infoframe.fnotes.bot.tnotes -relief sunken -bd 2 \
      -yscrollcommand "$infoframe.fnotes.bot.scroll set"
    scrollbar $infoframe.fnotes.bot.scroll -command "$infoframe.fnotes.bot tnotes yview"
    $infoframe.fnotes.bot.tnotes insert 1.0 $Subnet(Subnet${subnet_num}_notes)

    pack $infoframe.fnotes.top $infoframe.fnotes.bot -expand 0 -fill x -side top
    pack $infoframe.fnotes.top.lnotes -side left
    pack $infoframe.fnotes.bot -expand 1 -fill both -side top
    pack $infoframe.fnotes.bot.tnotes  -expand 1 -side left -fill both
    pack $infoframe.fnotes.bot.scroll -expand 0 -side left -fill y


    frame $w.fbuttons 
    button $w.fbuttons.ok -text "Done" -command "infoOk $w $subnet_num"
    button $w.fbuttons.clear -text "Clear All" \
      -command "infoClear $w $subnet_num"
    button $w.fbuttons.cancel -text "Cancel" \
      -command "infoCancel $w $subnet_num $backupDateTimeNotes"

    pack $infoframe.fname $infoframe.cdt $infoframe.fdt $infoframe.nv -side top -padx 1 -pady 1 -ipadx 2 -ipady 2 -fill x
    pack $infoframe.fnotes -expand 1 -fill both

    pack $w.fbuttons -side top -padx 1 -pady 1 -ipadx 2 -ipady 2 -fill x
    pack $w.fbuttons.ok $w.fbuttons.clear $w.fbuttons.cancel -side right -padx 5 -pady 5 -ipadx 3 -ipady 3
}

proc infoClear {w subnet_num} {
    global Subnet 
    set Subnet(Subnet${subnet_num}_userName) ""
    set Subnet(Subnet${subnet_num}_runDate) ""
    set Subnet(Subnet${subnet_num}_runTime) ""
    set Subnet(Subnet${subnet_num}_creationDate) ""
    set Subnet(Subnet${subnet_num}_creationTime) ""
    set Subnet(Subnet${subnet_num}_notes) ""
		set Subnet(Subnet${subnet_num}_relfilenames) ""
		set Subnet(Subnet${subnet_num}_netversion) ""
		
		set infoframe [$w.tab childsite 0]		
    $infoframe.fnotes.bot.tnotes delete 1.0 end
}

proc infoOk {w subnet_num} {
    global Subnet
		set infoframe [$w.tab childsite 0]
    set Subnet(Subnet${subnet_num}_notes) [$infoframe.fnotes.bot.tnotes get 1.0 end]
    networkHasChanged
    destroy $w
}

proc infoCancel {w subnet_num args } {
    global Subnet
    set Subnet(Subnet${subnet_num}_userName) [lindex $args 0]
    set Subnet(Subnet${subnet_num}_runDate) [lindex $args 1]
    set Subnet(Subnet${subnet_num}_runTime) [lindex $args 2]
    set Subnet(Subnet${subnet_num}_creationDate) [lindex $args 3]
    set Subnet(Subnet${subnet_num}_creationTime) [lindex $args 4]
    set Subnet(Subnet${subnet_num}_notes) [lindex $args 5]
		set Subnet(Subnet${subnet_num}_relfilenames) [lindex $args 6]
		set Subnet(Subnet${subnet_num}_netversion) [lindex $args 7]

    destroy $w
} 

proc loadfile {netedit_loadfile} {
    puts "NOTICE: `loadfile' has been disabled."
    puts "   To use old nets, remove the `loadfile' and `return' lines"
    puts "   from near the top of the file."
    return
}

proc loadnet { netedit_loadfile} {
    # Check to see of the file exists; warn user if it doesnt
    if { ![file exists $netedit_loadfile] } {
      set message "File \"$netedit_loadfile\" does not exist."
      createSciDialog -warning -message $message
      return
    }

    global netedit_savefile inserting PowerApp Subnet geometry
    if { !$inserting || ![string length $netedit_savefile] } {
      # Cut off the path from the net name and put in on the title bar:
      wm title . "SCIRun ([lindex [file split $netedit_loadfile] end])"
      # Remember the name of this net for future "Saves".
      set netedit_savefile $netedit_loadfile
    }

    renameSourceCommand
    
    # if we are not loading a powerapp network, show loading progress
    if { !$PowerApp } {
      showProgress 0 0 1 ;# -maybe- raise the progress meter
      setProgressText "Loading SCIRun network file..."
      update idletasks
      # The following counts the number of steps it will take to load the net
      resetScriptCount
      renameNetworkCommands counting_
      source $netedit_loadfile
      renameNetworkCommands counting_
      global scriptCount
      addProgressSteps $scriptCount(Total)
    }

    uplevel \#0 source \{$netedit_loadfile\}    
    resetSourceCommand

    if { !$PowerApp } {
      hideProgress
    }
    set Subnet(Subnet$Subnet(Loading)_Filename) $netedit_loadfile
    if { !$inserting } { setGlobal NetworkChanged 0 }
}

proc SCIRunNew_source { args } {
    set filename [lindex $args 0]
    if { [file exists $filename] } {
      return [uplevel 1 SCIRunBackup_source \{$filename\}]
    }

    set lastSettings [string last .settings $filename]
    if { ($lastSettings != -1) && \
	     ([expr [string length $filename] - $lastSettings] == 9) } {
      set file "[netedit getenv SCIRUN_SRCDIR]/nets/default.settings"
      global recentlyWarnedAboutDefaultSettings
      if { ![info exists recentlyWarnedAboutDefaultSettings] } {
          set recentlyWarnedAboutDefaultSettings 1
          after 10000 uplevel \#0 unset recentlyWarnedAboutDefaultSettings
          displayErrorWarningOrInfo "*** No $filename file was found.\n*** Loading the default dataset .settings file: $file\n" info
      }
      return [uplevel 1 SCIRunBackup_source \{$file\}]
    }
    puts "SCIRun TCL cannot source \'$filename\': File does not exist."
}

proc renameSourceCommand {} {
    if { [llength [info commands SCIRunBackup_source]] } return
    rename source SCIRunBackup_source
    rename SCIRunNew_source source
}



proc resetSourceCommand {} {
    if { ![llength [info commands SCIRunBackup_source]] } return
    rename source SCIRunNew_source
    rename SCIRunBackup_source source

}

proc showInvalidDatasetPrompt {} {
    set filename [lindex [file split [info script]] end]
    return [createSciDialog -warning \
		-parent . \
		-button1 "Select Dataset" \
		-button2 "Ignore Dataset" \
		-button3 "Exit SCIRun" \
		-message "This network can work on different datasets, but it\ncannot find files matching the pattern:\n\$(SCIRUN_DATA)/\$(SCIRUN_DATASET)/\$(SCIRUN_DATAFILE)*\n\nWhere:\nSCIRUN_DATA = [netedit getenv SCIRUN_DATA]\nSCIRUN_DATASET = [netedit getenv SCIRUN_DATASET]\nSCIRUN_DATAFILE = [netedit getenv SCIRUN_DATAFILE]\n\nThe environment variables listed are not set or they are invalid.\n\nChoose 'Select Dataset' to select a valid dataset directory.\neg: /usr/sci/data/SCIRunData/[netedit getenv SCIRUN_VERSION]/sphere\n\nChoose 'Ignore Dataset' to load $filename anyway.  You will\nhave to manually set the reader modules to valid filenames.\n\nIf you set these environment variables before running SCIRun, then\nyou will not receive this dialog when loading $filename.\n\nYou may also permanently set these variables in the ~/.scirunrc file."]
}



proc showChooseDatasetPrompt { initialdir } {
    if { ![string length $initialdir] } {
	set version [netedit getenv SCIRUN_VERSION]
	set initialdir "/usr/sci/data/SCIRunData/${version}"
	if { ![file exists $initialdir] } {
	    set initialdir [netedit getenv SCIRUN_OBJDIR]
	}
    }
    set value [tk_chooseDirectory -mustexist 1 -initialdir $initialdir \
		  -parent . -title "Select Dataset"]
    return $value
}


# sourceSettingsFile()
#
# Finds and then sources the ".settings" file.  Uses environment
# variables SCIRUN_DATA (for directory) and SCIRUN_DATASET (for data
# set) to find .settings file.  If these environment variables are not
# set or if they to not point to a valid file, then proc asks the user
# for input.  
#
# Returns "DATADIR DATASET"
#
proc sourceSettingsFile {} {
    # Attempt to get environment variables:
    set DATADIR [netedit getenv SCIRUN_DATA]
    set DATASET [netedit getenv SCIRUN_DATASET]
    set DATAFILE [netedit getenv SCIRUN_DATAFILE]
    
    if { ![string length $DATASET] } {
      # if env var SCIRUN_DATASET not set... default to sphere:
      set DATASET sphere
    } 

    global recentlyCalledSourceSettingsFile
    if { ![info exists recentlyCalledSourceSettingsFile] } {
      set initialdir ""
      while {![string length [glob -nocomplain "$DATADIR/$DATASET/$DATAFILE*"]] } {
          case [showInvalidDatasetPrompt] {
        1 "set data [showChooseDatasetPrompt $initialdir]"
        2 { 
            displayErrorWarningOrInfo "*** SCIRUN_DATA not set.  Reader modules will need to be manually set to valid filenames." warning
            uplevel \#0 source "[netedit getenv SCIRUN_SRCDIR]/nets/default.settings"
            return
        }
        3 "netedit quit"
	    }
	    if { [string length $data] } {
        set initialdir $data
        set data [file split $data]
        set DATASET [lindex $data end]
        set DATADIR [eval file join [lrange $data 0 end-1]]
	    }	      
	}

	displayErrorWarningOrInfo "*** Using SCIRUN_DATA=$DATADIR" info
	displayErrorWarningOrInfo "*** Using SCIRUN_DATASET=$DATASET" info
	displayErrorWarningOrInfo "*** Using SCIRUN_DATAFILE=$DATAFILE" info
	
	netedit setenv SCIRUN_DATA "$DATADIR"
	netedit setenv SCIRUN_DATASET "$DATASET"
	netedit setenv SCIRUN_DATAFILE "$DATAFILE"

	setGlobal recentlyCalledSourceSettingsFile 1
	after 10000 uplevel \#0 unset recentlyCalledSourceSettingsFile
    }

    set settings "$DATADIR/$DATASET/$DATASET.settings"
    uplevel 1 source $settings
    return "$DATADIR $DATASET $DATAFILE"
}

proc displayErrorWarningOrInfo { msg status } {
    AddToLog "$msg\n" "$status"
}


proc hideProgress { args } {
    if { ![winfo exists .splash] } return
    update idletasks
    update
    if {[llength $args]} {
      wm withdraw .splash
    } else {
      after 1500 { wm withdraw .splash }
    }
    grab release .splash
    update idletasks
}

proc showProgress { { show_image 0 } { steps none } { dismissButton 0 } { over . } } {

    global Color

    if { [envBool SCIRUN_HIDE_PROGRESS] && ![winfo exists .standalone] } return
    update idletasks
    set w .splash
    if { ![winfo exists $w] } {
      toplevel $w -bd 0 -relief raised -background $Color(MainBackGround)
      set width 820
      set height 100
      
      # do not show if either env vars are set to true, the only 
      # exception is when we are calling this from a powerapp's
      # show_help
      if { [envBool SCIRUN_NOSPLASH] || [envBool SCI_NOSPLASH] } {
        set show_image 0
      }      
      
      if { $show_image } {
        set height 600
      }
  
      
      set xoffset [expr ([winfo screenwidth $w] - $width)/2]
      set yoffset [expr ([winfo screenheight $w] - $height)/2]
      
      wm withdraw $w 
      wm geometry $w "${width}x${height}+${xoffset}+${yoffset}"
      wm protocol $w WM_DELETE_WINDOW hideProgress
    }
    if { $show_image } {
        wm title $w "Welcome to SCIRun v[netedit getenv SCIRUN_VERSION]"
    } else {
      wm title $w "Loading Network..."
    }

    if { ![winfo exists $w.frame] } {
      frame $w.frame -background $Color(MainBackGround)
      pack $w.frame -expand 1 -fill both
    }
    set w $w.frame


    if {[winfo exists .standalone]} {
      set show_image 1
    }

    if { [winfo exists $w.fb] } {
      destroy $w.fb
      unsetIfExists progressMeter
    }

    if { "$steps" != "none" } {
        # Hack: The more steps that are added to a feedback widget,
        # the more it (I assume) loses accuracy and thus causes the
        # display window to grow...  Therefore we are going to just use
        # 100 steps all the time.
      iwidgets::feedback $w.fb -steps 100 -barheight 11 -background $Color(MainBackGround)
      setGlobal progressMeter $w.fb
    }

    if { ![winfo exists $w.dismissBtn] } {
      button $w.dismissBtn -text "Dismiss" -width 10 \
        -command "hideProgress 1" -bd $Color(BorderWidth)\
        -activebackground $Color(MenuSelectBackGround)  \
        -activeforeground $Color(MenuSelectForeGround) \
        -background $Color(ButtonBackGround) \
        -foreground $Color(ButtonForeGround)  
    }

    if { $show_image } {
      global splashImageFile splashAboutFile
      if { [string length [info commands ::img::splash]] && \
               ![string equal $splashImageFile [::img::splash cget -file]] } {
          image delete ::img::splash
      }

      if { [string length [info commands ::img::about]] && \
               ![string equal $splashImageFile [::img::about cget -file]] } {
          image delete ::img::about
      }

      if { [file isfile $splashImageFile] && \
               [file readable $splashImageFile] && \
               ![string length [info commands ::img::splash]] } {
          image create photo ::img::splash -file $splashImageFile
      }

      if { [file isfile $splashAboutFile] && \
               [file readable $splashAboutFile] && \
               ![string length [info commands ::img::about]] } {
          image create photo ::img::about -file $splashAboutFile
      }
      if { ![winfo exists $w.splash] } {
          label $w.splash
      }
  
      if { $show_image == 2} {
        if { [string length [info command ::img::about]] } {
            $w.splash configure -image ::img::about -background $Color(MainBackGround)
            pack $w.splash -anchor w -fill both -expand yes
        } else {
            pack forget $w.splash
        }
      } else {
        if { [string length [info command ::img::splash]] } {
            $w.splash configure -image ::img::splash
            pack $w.splash -anchor w -fill both -expand yes
        } else {
            pack forget $w.splash
        }
      }
    } else {
      pack forget $w.splash
    }

    if { "$steps" != "none" } {
      setProgressText Initializing...
      addProgressSteps $steps
      pack $w.fb -padx 5 -fill x
    } else {
      pack forget $w.fb
    }

    if { $dismissButton } {
      pack $w.dismissBtn -side bottom -padx 5 -pady 3 -fill none
    } else {
      pack forget $w.dismissBtn
    }

    if {![winfo exists $w.scaffold] } {
      frame $w.scaffold -width 500 -background $Color(MainBackGround)
      pack $w.scaffold -side bottom
    }

    global PowerApp
    if { !$PowerApp } {
       wm deiconify .splash
       # centerWindow .splash $over
    } else {
        # PowerApp splash doesn't have a "." to center over...
        centerWindow .splash
    }
}

proc addProgressSteps { steps } {
    global totalProgressSteps
    set totalProgressSteps [expr $totalProgressSteps + $steps]
}

proc setProgressText { text } {
    global progressMeter
    if { [info exists progressMeter] && [winfo exists $progressMeter] } {
	$progressMeter configure -labeltext $text
	update idletasks
    }
}

proc incrProgress { { steps 1 } } {
    global progressMeter currentProgressStep currentPercent totalProgressSteps

    if { ![info exists progressMeter] } return

    set currentProgressStep [ expr $currentProgressStep + $steps ]
    # Covert (ie: the .0) to a floating point, and them multiply by 100 to get a number between 0-100.
    set percent             [expr (($currentProgressStep / $totalProgressSteps.0) * 100) - 1]

    ## The following debug statement is useful in figuring out the
    ## correct number of steps that the progress bar should be set to
    ## for the PowerApps.  When the number is determined, it should be
    ## updated in src/Dataflow/TCLThread/TCLThread.cc.
    #
    # puts "incrProgress $currentProgressStep"

    while { $currentPercent < $percent } {

        incr currentPercent

        if { $currentPercent >= 100 } {
            return
        }
	$progressMeter step
    }
}

    
	 

global LicenseResult
set licenseResult decline


proc acceptLicense {{ lic ""}} {

    netedit update_rcfile SCIRUN_ACCEPT_LICENSE ${lic}
}


proc licenseDialog { {firsttime 0} } {

    global licenseResult userData Color
    
    # Ignore license while regression testing
    if {[string compare [netedit getenv SCI_REGRESSION_TESTING] ""]} {
      return
    }
    
	if {[file exists [file join [netedit getenv SCIRUN_OBJDIR] Licenses]] } {
		set licensedir  [file join [netedit getenv SCIRUN_OBJDIR] Licenses]
	} else {
		if {[file exists [file join [netedit getenv SCIRUN_OBJDIR] .. Licenses]] } {
			set licensedir  [file join [netedit getenv SCIRUN_OBJDIR] .. Licenses]
		} else {
			# no license files available
			return
		}
	}
	
    if { $firsttime } {
      set lic ""
      foreach file [glob -nocomplain -directory $licensedir *.license] {
        set name [lindex [file split [file rootname $file]] end]
        set lic "${lic}_${name}"
      }

      set curlic [netedit getenv SCIRUN_ACCEPT_LICENSE]
      
      if {[string compare $curlic $lic] == 0} return
    }

    sci_toplevel .license
    wm title .license {SCIRun License Information}
    sci_tabnotebook .license.tabs -height 400 -width 600 -tabpos n

    .license.tabs add -label "Licenses"
    set firsttab [.license.tabs childsite 0]

    sci_labeledframe $firsttab.info -labeltext "License Information"
      
    set frame [$firsttab.info childsite]
    
    sci_scrolledtext $frame.text -vscrollmode dynamic -hscrollmode dynamic -visibleitems 70x30
    set textw [$frame.text component text]
    pack $firsttab.info -fill both -expand 1
    pack $frame.text -fill both -expand 1
    
    set sr_version [netedit getenv SCIRUN_VERSION]
    $textw insert end "SCIRun $sr_version License Information\n\n"
    $textw insert end "This version of SCIRun includes the following external packages:\n\n"
    
    .license.tabs select 0

    set i 1

    set file "$licensedir/SCIRun.license"
    set name [lindex [file split [file rootname $file]] end]
    
    .license.tabs add -label "$name License"
    set tab [.license.tabs childsite $i]

    sci_labeledframe $tab.text -labeltext "$name License"
    pack $tab.text -fill both -expand 1

    set text [$tab.text childsite]
    sci_scrolledtext $text.license -vscrollmode dynamic -hscrollmode dynamic -visibleitems 70x30
    
    $text.license import $file
    pack $text.license -fill both -expand 1

    foreach file [glob -nocomplain -directory $licensedir *.license] {
      set name [lindex [file split [file rootname $file]] end]
              
      if {[string compare $name "SCIRun"]} {
        incr i
        .license.tabs add -label "$name License"
        set tab [.license.tabs childsite $i]

        sci_labeledframe $tab.text -labeltext "$name License"
        pack $tab.text -fill both -expand 1

        set text [$tab.text childsite]
        sci_scrolledtext $text.license -vscrollmode dynamic -hscrollmode dynamic -visibleitems 70x30
      
        $text.license import $file
        pack $text.license -fill both -expand 1

        $textw insert end " - ${name}\n"
      }
      
    }
    


    $textw insert end "\n\nOther SCIRun Binary and Source distributions are available at:\nhttp://software.sci.utah.edu"
    if { $firsttime } {
      $textw insert end "\n\nPlease read the license for each package before continuing.\n"
      $textw insert end "and press 'Accept licenses' to continue or 'Decline' to quit\n"
    } else {
      $textw insert end "\n\nThe license of SCIRun and the additional packages can be found by\nclicking on the different license tabs\n"
    }
     
    pack .license.tabs -side top -fill both -expand yes

    frame .license.buttons -background $Color(MainBackGround) -bd 1 -relief raised
    pack .license.buttons -side bottom -fill x 

    if { $firsttime } {
      set w .license.buttons
      sci_button $w.accept -text "Accept Licenses" -width 20 -command "netedit update_env SCIRUN_ACCEPT_LICENSE $lic; acceptLicense $lic; catch {destroy .license}"
      sci_button $w.decline -text "Decline" -width 12 \
          -command "netedit quit"
      pack $w.accept $w.decline   -padx 5 -pady 5 -side right -padx 20

      wm protocol .license WM_DELETE_WINDOW {
          createSciDialog -error -message "You must choose Accept or Decline to continue."}
    } else {
      wm protocol .license WM_DELETE_WINDOW { destroy .license }
      
      set w .license.buttons
      sci_button $w.ok -text "OK" -width 10 -command {destroy .license} 
      pack $w.ok -padx 2 -pady 2 -side bottom
    
    }

    moveToCursor .license
    wm deiconify .license
    grab .license
    if { $firsttime } { tkwait window .license }    
    return

}


proc promptUserToCopySCIRunrc {} {
    global copyResult Color
    set w .copyRCprompt

    sci_toplevel $w
    wm withdraw $w

    set copyResult 0
    set dontAskAgain 0

    set sr_version [netedit getenv SCIRUN_VERSION]
    set rc_version [netedit getenv SCIRUN_RCFILE_VERSION]

    if { $rc_version == "" } {
        set rc_version "bak"
    }

    wm title $w "Update .scirunrc file to version $sr_version?"
    sci_frame $w.ff
    sci_label $w.ff.message \
     -text "A newer version of the .scirunrc file is avaliable with this release.\n\nThis file contains SCIRun environment variables that are necessary\nfor new features.\n\nPlease note: If you have made changes to your .scirunrc file they will\nbe undone by updating to a new .scirunrc.\n\nHowever, if you generate a new .scirunrc, your existing file will be\nbacked up to .scirunrc.$rc_version just in case.\n\nWould you like SCIRun to update to the new .scirunrc?" -justify left

    frame $w.but -background $Color(MainBackGround)
    sci_button $w.but.ok -text "Update" -command "set copyResult 1" 
     
    sci_button $w.but.no -text "Don't Update" -command "set copyResult 0"
    sci_button $w.but.dontask -text "Don't Ask Again" -command "set copyResult 2"  


    pack $w.but.ok $w.but.no $w.but.dontask  -side left -pady 3 -padx 5 -expand 1
    pack $w.ff -expand 1 -fill both
    pack $w.ff.message -expand 1 -fill both -padx 10 -pady 10

    # Vertical separator
    # frame $w.separator -height 2 -relief sunken -borderwidth 2 -background $Color(MainBackGround)
    # pack  $w.separator -fill both -padx 5 -pady 1 -expand no

    pack $w.but -expand 1 -fill both -side top -anchor n

    # Override the destroy window decoration and make it not do anything
    wm protocol $w WM_DELETE_WINDOW "SciRaise $w"
    moveToCursor $w
    SciRaise $w

    vwait copyResult

    if { $copyResult == 2 } {
      if { [catch { set rcfile [open ~/.scirunrc "WRONLY APPEND"] }] } {
                puts "Unable to open ~/.scirunrc"
                return 0
            }

      puts $rcfile \
                "\n\# This section added when the user chose 'Dont Ask Again'"
      puts $rcfile \
                "\# when prompted about updating the .scirurc file version"

      set rc_sub_version [netedit getenv SCIRUN_RCFILE_SUBVERSION]
      puts $rcfile "SCIRUN_RCFILE_VERSION=${sr_version}.${rc_sub_version}"
      close $rcfile
      set copyResult 0
    }

    destroy $w
    return $copyResult
}
    
proc validFile { args } {
    set name [lindex $args 0]
    return [expr [file isfile $name] && [file readable $name]]
}

proc validDir { name } {
    return [expr [file isdirectory $name] && \
		 [file writable $name] && [file readable $name]]
}


# Removes the element at pos from a list without a set - similar to lappend
# ex: 
#   set bob "0 1 2 3 4 5"
#   listRemove bob 2
#   puts $bob
# output:
#   0 1 3 4 5
proc listRemove { name pos } {
    uplevel 1 set $name \[list [lreplace [uplevel 1 set $name] $pos $pos]\]
}

# Finds the first instance of elem in a list then removes it from the list 
# ex: 
#   set bob "foo bar foo2 bar2 foo3 bar3"
#   listFindAndRemove bob foo2
#   puts $bob
# output:
#   foo bar bar2 foo3 bar3
proc listFindAndRemove { name elem } {
    set elements [uplevel 1 set $name]
    set pos [lsearch $elements $elem]    
    uplevel 1 set $name \[list [lreplace $elements $pos $pos]\]
}


proc initVarStates { var save substitute } {
    set var [string trimleft $var :]
    if { [string first msg_stream $var] != -1 } return
    global ModuleSavedVars ModleSubstitutedVars
    set ids [split $var -]
    if { [llength $ids] < 2 } return
    set module [lindex $ids 0]
    set varname [join [lrange $ids 1 end] -]
    
    if { $save } {
      lappend ModuleSavedVars($module) $varname
      # TODO: find a faster mechanism than the next line for setting changed
      uplevel \#0 trace variable \"$var\" w networkHasChanged
    }
     
    if { $substitute } {
      lappend ModuleSubstitutedVars($module) $varname
    }
}

proc setVarStates { var save substitute isfilename} {

    global ModuleSavedVars ModuleSubstitutedVars ModuleIsFilenameVars
    set var [string trimleft $var :]
    if { [string first msg_stream $var] != -1 } return
    set ids [split $var -]
    set module [lindex $ids 0]
    set varname [join [lrange $ids 1 end] -]
    
    if { [expr [string length $varname] == 0] || [expr [string length $module] == 0 ] } {
      return
    }
    
    if { ![info exists ModuleSavedVars($module)] } {
			set saved 0
    } else {
			set saved [expr [lsearch $ModuleSavedVars($module) $varname] != -1]
    }

    if { $save && !$saved} {
			lappend ModuleSavedVars($module) $varname
			uplevel \#0 trace variable \"$var\" w networkHasChanged
    } elseif { !$save && $saved } {
			listFindAndRemove ModuleSavedVars($module) $varname
			uplevel \#0 trace vdelete \"$var\" w networkHasChanged
    }
          	
    if { ![info exists ModuleSubstitutedVars($module)] } {
			set substituted 0
    } else {
			set substituted [expr [lsearch $ModuleSubstitutedVars($module) $varname] != -1]
    }

    if { $substitute && !$substituted } {
			lappend ModuleSubstitutedVars($module) $varname
    } elseif { !$substitute && $substituted } {
			listFindAndRemove ModuleSubstitutedVars($module) $varname
    }

    if { ![info exists ModuleIsFilenameVars($module)] } {
			set filename 0
    } else {
			set filename [expr [lsearch $ModuleIsFilenameVars($module) $varname] != -1]
    }

    if { $isfilename && !$filename } {
			lappend ModuleIsFilenameVars($module) $varname
    } elseif { !$isfilename && $filename } {
			listFindAndRemove ModuleIsFilenameVars($module) $varname
    }
}


# Debug procedure to print global variable values
proc printvars { pattern } {
    foreach name [lsort [uplevel \#0 "info vars *${pattern}*"]] { 
      upvar \#0 $name var
      if { [uplevel \#0 array exists \"$name\"] } {
          uplevel \#0 parray \"$name\"
      } else {
          puts "set \"$name\" \{$var\}"
      }
    }
}

proc setGlobal { var val } {
    uplevel \#0 set \"$var\" \{$val\}
}

# Will only set global variable $var to $val if it doesnt already exist
proc initGlobal { var val } {
    upvar \#0 $var globalvar
    if { ![info exists globalvar] } {
      set globalvar $val
    }
}


proc popFront { listname } {
    upvar 1 $listname list
    if ![info exists list] return
    set front [lindex $list 0]
    set list [lrange $list 1 end]
    return $front
}

proc popBack { listname } {
    upvar 1 $listname list
    if ![info exists list] return
    set back [lindex $list end]
    set list [lrange $list 0 end-1]
    return $back
}

proc maybeWrite_init_DATADIR_and_DATASET { out } {
    global ModuleSubstitutedVars
    if { ![envBool SCIRUN_NET_SUBSTITUTE_DATADIR] } return
    foreach module [array names ModuleSubstitutedVars] {
      foreach var $ModuleSubstitutedVars($module) {
          upvar \#0 $module-$var val
          if { [info exists val] && \
        ![string equal $val [subDATADIRandDATASET $val]] } {
        netedit net-add-env-var scisub_datadir SCIRUN_DATA
        netedit net-add-env-var scisub_datafile SCIRUN_DATAFILE
        netedit net-add-env-var scisub_dataset SCIRUN_DATASET
        return
          }
      }
    }
}


proc maybeWriteTCLStyleCopyright { out } {
    if { ![envBool SCIRUN_INSERT_NET_COPYRIGHT] } return 
    catch "set license [open [netedit getenv SCIRUN_SRCDIR]/LICENSE]"
    if { ![info exists license] } return
    while { ![eof $license] } {
	puts $out "\# [gets $license]"
    }
    close $license
}


proc init_DATADIR_and_DATASET {} {
    uplevel 1 sourceSettingsFile
    upvar 1 DATADIR datadir DATASET dataset DATAFILE datafile
    set datadir [netedit getenv SCIRUN_DATA]
    set dataset [netedit getenv SCIRUN_DATASET]
    set datafile [netedit getenv SCIRUN_DATAFILE]
    netedit setenv SCIRUN_NET_SUBSTITUTE_DATADIR true
}
    
proc backupNetwork { } {
  # check if we have a filename and if so save it to #filename#
  global netedit_savefile
  if { [file exists $netedit_savefile] } {
      set root_filename [lindex [file split $netedit_savefile] end]

      # don't save back ups of back up files as ##filename##!!!
      if {[string index $root_filename 0] == "#" && 
	  [string index $root_filename end] == "#"} {
	      set root_filename [string range $root_filename 1 end-1]
      }

      set current_path [lrange [file split $netedit_savefile] 0 end-1]
      if {[llength $current_path] > 0} {
	  set current_path [eval file join [lrange [file split $netedit_savefile] 0 end-1]]
      } else {
	  # net in current directory so we need something to the path
	  # so the file joins work
	  set current_path [pwd]
      }

      # First attempt to save it to same directory as netedit_savefile
      if {[file writable $current_path]} {
	  set dest [eval file join $current_path \#$root_filename\#]
	  writeNetwork $dest
      } else {
	  # else save to home/SCIRun directory as nedetit_savefile
	  set src "[file split [netedit getenv HOME]]"
	  set src [lappend src SCIRun \#$root_filename\#]
	  set dest [eval file join $src]
	  writeNetwork $dest
      }
  } else {
      # Attempt to write to #MyNetwork.srn# in current directory
      set src  "[file split [pwd]] MyNetwork.srn"
      set src [eval file join $src]
      set parts [file split $src]
      set parts [lreplace $parts end end \#[lindex $parts end]\#]
      set dir [pwd]

      if {[file writable $dir]} {
          set dest [eval file join $parts]
	  writeNetwork $dest
      } else {
	  # Else write to home/SCIRun directory
	  set src "[file split [netedit getenv HOME]]"
	  set src [lappend src SCIRun \#MyNetwork.srn\#]
	  set dest [eval file join $src]
	  writeNetwork $dest
      }
  }
}

proc writeNetwork { filename { subnet 0 } } {
    # if the file already exists, back it up to "filename.bak"
    if { [file exists $filename] &&
         [string index $filename 0] != "#" &&
         [string index $filename end] != "#"} {
      set src  "[file split [pwd]] [file split ${filename}]"
      set src [eval file join $src]
      set parts [file split $src]
      set parts "[lreplace $parts end end [lindex $parts end]].bak"
      set dest [eval file join $parts]
      catch [file rename -force $src $dest]
    }

    global Subnet
    set netversion [set Subnet(Subnet0_netversion)]
    incr netversion
    set Subnet(Subnet0_netversion) $netversion
    netedit start-net-doc $filename v[netedit getenv SCIRUN_VERSION] \
      $netversion
    set out stdout
    
    #maybeWriteTCLStyleCopyright $out
    maybeWrite_init_DATADIR_and_DATASET $out
    genSubnetScript $subnet ""

    netedit write-net-doc
}


# Numerically compares two version strings and returns:
#  -1 if ver1 < ver2,
#   0 if ver1 == ver2
#   1 if ver1 > ver2
proc compareVersions { ver1 ver2 } {
    set v1 [split $ver1 .]
    set v2 [split $ver2 .]
    set l1 [llength $v1]
    set l2 [llength $v2]
    set len [expr ($l1 > $l2) ? $l1 : $l2]
    for { set i 0 } {$i < $len} {incr i} {
	set n1 -1
	set n2 -1
	if {$i < $l1} {
	    set n1 [lindex $v1 $i]
	}
	if {$i < $l2} {
	    set n2 [lindex $v2 $i]
	}
	if { $n1 < $n2 } {
	    return -1
	}
	if { $n2 < $n1 } {
	    return 1
	}
    }
    return 0
}

proc txt { args } {
    return [join $args \n]
}

proc in_power_app {} {
    return [winfo exists .standalone]
}

proc editrcfile {} {

  global Color
  set rcfile [netedit getenv "SCIRUN_RCFILE"]

  # Create a unique name for the file selection window
  set w .ui-edit-scirunrc

  # if the file selector is open, bring it to the front
  # in case it is iconified, deiconify
  if { [winfo exists $w] } {
      SciRaise $w
      return
  }

  toplevel $w -class TkFDialog        

  iwidgets::scrolledtext $w.file -labeltext "SCIRun Configuration file" \
    -visibleitems 70x20 -vscrollmode dynamic -hscrollmode dynamic -wrap none \
    -background $Color(MainBackGround)
 
  pack $w.file -fill both -expand true
  
  set b [set rcfile]
  $w.file import [set rcfile]
  frame $w.f1 -background $Color(MainBackGround)
  pack $w.f1 -side bottom -anchor e -fill x

  button $w.f1.load -text "Load" -command "$w.file clear; $w.file import $b" \
     -bd $Color(BorderWidth) \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(ButtonBackGround) \
     -foreground $Color(ButtonForeGround)    
  button $w.f1.save -text "Save" -command "$w.file export $b; netedit reload_rcfile" \
     -bd $Color(BorderWidth) \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(ButtonBackGround) \
     -foreground $Color(ButtonForeGround)   
  button $w.f1.done -text "Done" -command "destroy $w" \
     -bd $Color(BorderWidth) \
     -activebackground $Color(MenuSelectBackGround)  \
     -activeforeground $Color(MenuSelectForeGround) \
     -background $Color(ButtonBackGround) \
     -foreground $Color(ButtonForeGround)     
  
  frame $w.f1.space -background $Color(MainBackGround) -width 20
  pack $w.f1.space $w.f1.done $w.f1.save $w.f1.load  -anchor e -side right -padx 3 -pady 1
                    
}

# Some how info exists fails in the newest TCL version
# this function works around that bug
proc info_exists { varname } {
  upvar 1 $varname testvar
  if { [catch { set test [set testvar] } ] } {
    return 0
  } else {
    return 1          
  }
}

