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


set tcl_prompt1 "puts -nonewline \"scirun> \""
set tcl_prompt2 "puts -nonewline \"scirun>> \""

global port_spacing
set port_spacing 14

global port_width
set port_width 11

global port_height
set port_height 5 

global port_light_height
set port_light_height 4

global smod_button_height
set smod_button_height 16

global lmod_button_height
set lmod_button_height 22


global Color
set Color(Selected) LightSkyBlue2
set Color(Disabled) black
set Color(Compiling) "\#f0e68c"
set Color(Trace) red
set Color(ConnDisabled) gray
set Color(NetworkEditor) "\#161628"
set Color(NetworkEditorSmall) "\#101028"
set Color(SubnetEditor) "\#250030"
set Color(ErrorFrameBG) $Color(NetworkEditor)
set Color(ErrorFrameFG) white
set Color(IconFadeStart) $Color(NetworkEditor)
set Color(Basecolor) "\#cccccc"
set Color(MenuSelectBackGround) "\#808090"
set Color(MenuSelectForeGround) white
set Color(MenuSelectBackGroundHL) "\#ff1111"
set Color(MenuSelectForeGroundHL) "\#e5e5cc"
set Color(MenuBackGround) "\#eeeef3"
set Color(MenuForeGround) black
set Color(MenuHeadingForeGround) "\#000099"
set Color(MenuHeadingBackGround) "\#eeeef3"
set Color(MenuBackGroundHL) "\#e5e5cc"
set Color(MenuForeGroundHL) "\#ff1111"
set Color(MenuBarBackGround) "\#bbbbd4"
set Color(MenuBarForeGround) black
set Color(Trough) "\#aaaaaa"
set Color(ButtonBackGround) "\#d4d4dd"
set Color(ButtonForeGround) black
set Color(MainBackGround) "\#afafb4"
set Color(MainForeGround) black
set Color(EditBackGround) "\#ffffff"
set Color(EditHighLight) "\#ffaaaa"
set Color(EditForeGround) black
set Color(UIBackDrop) "\#d8d8de"
set Color(UIBackGround) "\#e5e5e9"
set Color(UIForeGround) black

set Color(LockColor) "\#aa0000"

set Color(SCIText) "\#444444"
set Color(Module) "\#cccccc"
set Color(ModuleButton) "\#bbbbbb"
set Color(ModuleError) "\#ee1111"
set Color(ModuleWarning) yellow
set Color(ModuleRemark) "\#4444ee"
set Color(ModuleProgress) "\#00aa00"

set Color(Waiting) LightGoldenrod3
set Color(Executing) "\#aaccaa"
set Color(BorderWidth) 1


global Subnet
set Subnet(Loading) 0

set basecolor $Color(Basecolor)

. configure -background $basecolor

option add *Frame*background black

option add *Button*padX 1
option add *Button*padY 1

option add *background $basecolor
option add *activeBackground $basecolor
option add *sliderForeground $basecolor
option add *troughColor $basecolor
option add *activeForeground white

option add *Scrollbar*activeBackground $basecolor
option add *Scrollbar*foreground $basecolor
option add *Scrollbar*width .35c
option add *Scale*width .35c

option add *selectBackground "white"
option add *selectForeground "black"
option add *selector red
option add *font "-Adobe-Helvetica-bold-R-Normal--*-120-75-*"
option add *highlightThickness 0

