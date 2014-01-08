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


itcl::class SCIRun_Visualization_CreateStandardColorMaps { 
    inherit Module 
    protected variable exposed
    protected variable colorMaps
    protected variable colorMap
    protected variable alphaMap
    protected variable curX
    protected variable curY
    protected variable selected
     constructor { {args ""} } {
        eval configure $args 
        set name CreateStandardColorMaps 
        set_defaults 
	buildColorMaps
    } 
    
    method set_defaults {} { 
	set exposed 0
	set colorMap {}
	set selected -1

	trace variable $this-mapType w "$this lookupOldIndex"
    }   
    
    method buildColorMaps {} {
	set colorMaps {
	    { "Gray" { { 0 0 0 } { 51 51 51 } { 102 102 102 } { 153 153 153 } { 204 204 204 } { 255 255 255 } } }
	    { "Rainbow" {
		{0 0 255} {0 52 255}
		{1 80 255} {3 105 255}
		{5 132 255} {9 157 243}
		{11 177 213} {15 193 182}
		{21 210 152} {30 225 126}
		{42 237 102} {60 248 82}
		{87 255 62} {116 255 49}
		{148 252 37} {178 243 27}
		{201 233 19} {220 220 14}
		{236 206 10} {247 185 8}
		{253 171 5} {255 151 3}
		{255 130 2} {255 112 1}
		{255 94 0} {255 76 0}
		{255 55 0} {255 0 0}}}
	    { "Old Rainbow" {
		{ 0 0 255}   { 0 102 255}
		{ 0 204 255}  { 0 255 204}
		{ 0 255 102}  { 0 255 0}
		{ 102 255 0}  { 204 255 0}
		{ 255 234 0}  { 255 204 0}
		{ 255 102 0}  { 255 0 0}}}
	    { "Darkhue" {
		{ 0  0  0 }  { 0 28 39 }
		{ 0 30 55 }  { 0 15 74 }
		{ 1  0 76 }  { 28  0 84 }
		{ 32  0 85 }  { 57  1 92 }
		{ 108  0 114 }  { 135  0 105 }
		{ 158  1 72 }  { 177  1 39 }
		{ 220  10 10 }  { 229 30  1 }
		{ 246 72  1 }  { 255 175 36 }
		{ 255 231 68 }  { 251 255 121 }
		{ 239 253 174 }}}
	    { "Lighthue" {
		{ 64  64  64 }  { 64 80 84 }
		{ 64 79 92 }  { 64 72 111 }
		{ 64  64 102 }  { 80 64 108 }
		{ 80 64 108 }  { 92  64 110 }
		{ 118  64 121 }  { 131  64 116 }
		{ 133  64 100 }  { 152  64 84 }
		{ 174  69 69 }  { 179 79  64 }
		{ 189 100  64 }  { 192 152 82 }
		{ 192 179 98 }  { 189 192 124 }
		{ 184 191 151 }}}
	    { "Blackbody" {
		{0 0 0}   {52 0 0}
		{102 2 0}   {153 18 0}
		{200 41 0}   {230 71 0}
		{255 120 0}   {255 163 20}
		{255 204 55}   {255 228 80}
		{255 247 120}   {255 255 180}
		{255 255 255}}}
	    { "Don" {
		{   0  90 255 }    {  51 104 255 }
		{ 103 117 255 }    { 166 131 245 }
		{ 181 130 216 }    { 192 129 186 }
		{ 197 128 172 }    { 230 126  98 }
		{ 240 126  49 }    { 255 133   0 }}}
	    { "BP Seismic" { { 0 0 255 } { 255 255 255} { 255 0 0 } } }
	    { "Dark Gray" {
		{   0  0  0 }    {  0 0 0 }
		{ 128 128 128 } { 255 255 255 }}}
	    { "Red Tint" { { 20 0 0 } { 255 235 235 } } }
	    { "Orange Tint" { { 20 10 0 } { 255 245 235 } } }
	    { "Yellow Tint" { { 20 20 0 } { 255 255 235 } } }
	    { "Green Tint" { { 0 20 0 } { 235 255 235 } } }
	    { "Cyan Tint" { { 0 20 20 } { 235 255 255 } } }
	    { "Blue Tint" { { 0 0 20 } { 235 235 255 } } }
	    { "Purple Tint" { { 10 0 20 } { 245 235 255 } } }
	    { "Pink,White,Blue" { 
		{ 204  0  153 }   { 206  15 159 }
		{ 208  32 164 }   { 210  48 170 }
		{ 211  65 175 }   { 213  81 181 }
		{ 215  98 186 }   { 217 114 192 }
		{ 219 131 197 }   { 221 147 203 }
		{ 223 164 208 }   { 224 180 214 }
		{ 226 196 219 }   { 228 213 225 }
                { 255 255 255 }   { 213 213 255 }
		{ 196 196 255 }   { 180 180 255 }
		{ 164 164 255 }   { 147 147 255 }
		{ 131 131 255 }   { 114 114 255 }
		{ 98  98  255 }   { 81  81  255 }
		{ 65  65  255 }   { 48  48  255 }
		{ 32  32  255 }   { 15  15  255 }
		{ 0   0   255 } } }
	    { "Orange,Black,Lime" { 
		{ 255 153   0 }   { 240 145   0 } 
		{ 224 137   0 }   { 208 128   0 } 
		{ 192 119   0 }   { 176 111   0 } 
		{ 160 102   0 }   { 144  93   0 } 
		{ 129  85   0 }   { 113  76   0 } 
		{  97  68   0 }   {  81  59   0 } 
		{  65  50   0 }   {  49  42   0 } 
                {  32  32  32 }   {  42  49   0 } 
		{  50  65   0 }   {  59  81   0 } 
		{  68  97   0 }   {  76 113   0 } 
		{  85 129   0 }   {  93 144   0 } 
		{ 102 160   0 }   { 111 176   0 } 
		{ 119 192   0 }   { 128 208   0 }   
		{ 137 224   0 }   { 145 240   0 }   
      { 153 255   0 } } }
{ "Mixed Rainbow" {
{0 0 131}
{0 68 255}
{4 255 251}
{195 255 60}
{255 124 0}
{187 0 0}
{0 0 135}
{0 72 255}
{8 255 247}
{199 255 56}
{255 120 0}
{183 0 0}
{0 0 139}
{0 76 255}
{12 255 243}
{203 255 52}
{255 116 0}
{179 0 0}
{0 0 143}
{0 80 255}
{16 255 239}
{207 255 48}
{255 112 0}
{175 0 0}
{0 0 147}
{0 84 255}
{20 255 235}
{211 255 44}
{255 108 0}
{171 0 0}
{0 0 151}
{0 88 255}
{24 255 231}
{215 255 40}
{255 104 0}
{167 0 0}
{0 0 155}
{0 92 255}
{28 255 227}
{219 255 36}
{255 100 0}
{163 0 0}
{0 0 159}
{0 96 255}
{32 255 223}
{223 255 32}
{255 96 0}
{159 0 0}
{0 0 163}
{0 100 255}
{36 255 219}
{227 255 28}
{255 92 0}
{155 0 0}
{0 0 167}
{0 104 255}
{40 255 215}
{231 255 24}
{255 88 0}
{151 0 0}
{0 0 171}
{0 108 255}
{44 255 211}
{235 255 20}
{255 84 0}
{147 0 0}
{0 0 175}
{0 112 255}
{48 255 207}
{239 255 16}
{255 80 0}
{143 0 0}
{0 0 179}
{0 116 255}
{52 255 203}
{243 255 12}
{255 76 0}
{139 0 0}
{0 0 183}
{0 120 255}
{56 255 199}
{247 255 8}
{255 72 0}
{135 0 0}
{0 0 187}
{0 124 255}
{60 255 195}
{251 255 4}
{255 68 0}
{131 0 0}
{0 0 191}
{0 128 255}
{64 255 191}
{255 255 0}
{255 64 0}
{128 0 0}
{0 0 195}
{0 131 255}
{68 255 187}
{255 251 0}
{255 60 0}
{0 0 199}
{0 135 255}
{72 255 183}
{255 247 0}
{255 56 0}
{0 0 203}
{0 139 255}
{76 255 179}
{255 243 0}
{255 52 0}
{0 0 207}
{0 143 255}
{80 255 175}
{255 239 0}
{255 48 0}
{0 0 211}
{0 147 255}
{84 255 171}
{255 235 0}
{255 44 0}
{0 0 215}
{0 151 255}
{88 255 167}
{255 231 0}
{255 40 0}
{0 0 219}
{0 155 255}
{92 255 163}
{255 227 0}
{255 36 0}
{0 0 223}
{0 159 255}
{96 255 159}
{255 223 0}
{255 32 0}
{0 0 227}
{0 163 255}
{100 255 155}
{255 219 0}
{255 28 0}
{0 0 231}
{0 167 255}
{104 255 151}
{255 215 0}
{255 24 0}
{0 0 235}
{0 171 255}
{108 255 147}
{255 211 0}
{255 20 0}
{0 0 239}
{0 175 255}
{112 255 143}
{255 207 0}
{255 16 0}
{0 0 243}
{0 179 255}
{116 255 139}
{255 203 0}
{255 12 0}
{0 0 247}
{0 183 255}
{120 255 135}
{255 199 0}
{255 8 0}
{0 0 251}
{0 187 255}
{124 255 131}
{255 195 0}
{255 4 0}
{0 0 255}
{0 191 255}
{128 255 128}
{255 191 0}
{255 0 0}
{0 4 255}
{0 195 255}
{131 255 124}
{255 187 0}
{251 0 0}
{0 8 255}
{0 199 255}
{135 255 120}
{255 183 0}
{247 0 0}
{0 12 255}
{0 203 255}
{139 255 116}
{255 179 0}
{243 0 0}
{0 16 255}
{0 207 255}
{143 255 112}
{255 175 0}
{239 0 0}
{0 20 255}
{0 211 255}
{147 255 108}
{255 171 0}
{235 0 0}
{0 24 255}
{0 215 255}
{151 255 104}
{255 167 0}
{231 0 0}
{0 28 255}
{0 219 255}
{155 255 100}
{255 163 0}
{227 0 0}
{0 32 255}
{0 223 255}
{159 255 96}
{255 159 0}
{223 0 0}
{0 36 255}
{0 227 255}
{163 255 92}
{255 155 0}
{219 0 0}
{0 40 255}
{0 231 255}
{167 255 88}
{255 151 0}
{215 0 0}
{0 44 255}
{0 235 255}
{171 255 84}
{255 147 0}
{211 0 0}
{0 48 255}
{0 239 255}
{175 255 80}
{255 143 0}
{207 0 0}
{0 52 255}
{0 243 255}
{179 255 76}
{255 139 0}
{203 0 0}
{0 56 255}
{0 247 255}
{183 255 72}
{255 135 0}
{199 0 0}
{0 60 255}
{0 251 255}
{187 255 68}
{255 131 0}
{195 0 0}
{0 64 255}
{0 255 255}
{191 255 64}
{255 128 0}
{191 0 0}}}
{ "Mixed GrayScale" {
{0 0 0}
{48 48 48}
{96 96 96}
{144 144 144}
{192 192 192}
{240 240 240}
{1 1 1}
{49 49 49}
{97 97 97}
{145 145 145}
{193 193 193}
{241 241 241}
{2 2 2}
{50 50 50}
{98 98 98}
{146 146 146}
{194 194 194}
{242 242 242}
{3 3 3}
{51 51 51}
{99 99 99}
{147 147 147}
{195 195 195}
{243 243 243}
{4 4 4}
{52 52 52}
{100 100 100}
{148 148 148}
{196 196 196}
{244 244 244}
{5 5 5}
{53 53 53}
{101 101 101}
{149 149 149}
{197 197 197}
{245 245 245}
{6 6 6}
{54 54 54}
{102 102 102}
{150 150 150}
{198 198 198}
{246 246 246}
{7 7 7}
{55 55 55}
{103 103 103}
{151 151 151}
{199 199 199}
{247 247 247}
{8 8 8}
{56 56 56}
{104 104 104}
{152 152 152}
{200 200 200}
{248 248 248}
{9 9 9}
{57 57 57}
{105 105 105}
{153 153 153}
{201 201 201}
{249 249 249}
{10 10 10}
{58 58 58}
{106 106 106}
{154 154 154}
{202 202 202}
{250 250 250}
{11 11 11}
{59 59 59}
{107 107 107}
{155 155 155}
{203 203 203}
{251 251 251}
{12 12 12}
{60 60 60}
{108 108 108}
{156 156 156}
{204 204 204}
{252 252 252}
{13 13 13}
{61 61 61}
{109 109 109}
{157 157 157}
{205 205 205}
{253 253 253}
{14 14 14}
{62 62 62}
{110 110 110}
{158 158 158}
{206 206 206}
{254 254 254}
{15 15 15}
{63 63 63}
{111 111 111}
{159 159 159}
{207 207 207}
{255 255 255}
{16 16 16}
{64 64 64}
{112 112 112}
{160 160 160}
{208 208 208}
{17 17 17}
{65 65 65}
{113 113 113}
{161 161 161}
{209 209 209}
{18 18 18}
{66 66 66}
{114 114 114}
{162 162 162}
{210 210 210}
{19 19 19}
{67 67 67}
{115 115 115}
{163 163 163}
{211 211 211}
{20 20 20}
{68 68 68}
{116 116 116}
{164 164 164}
{212 212 212}
{21 21 21}
{69 69 69}
{117 117 117}
{165 165 165}
{213 213 213}
{22 22 22}
{70 70 70}
{118 118 118}
{166 166 166}
{214 214 214}
{23 23 23}
{71 71 71}
{119 119 119}
{167 167 167}
{215 215 215}
{24 24 24}
{72 72 72}
{120 120 120}
{168 168 168}
{216 216 216}
{25 25 25}
{73 73 73}
{121 121 121}
{169 169 169}
{217 217 217}
{26 26 26}
{74 74 74}
{122 122 122}
{170 170 170}
{218 218 218}
{27 27 27}
{75 75 75}
{123 123 123}
{171 171 171}
{219 219 219}
{28 28 28}
{76 76 76}
{124 124 124}
{172 172 172}
{220 220 220}
{29 29 29}
{77 77 77}
{125 125 125}
{173 173 173}
{221 221 221}
{30 30 30}
{78 78 78}
{126 126 126}
{174 174 174}
{222 222 222}
{31 31 31}
{79 79 79}
{127 127 127}
{175 175 175}
{223 223 223}
{32 32 32}
{80 80 80}
{128 128 128}
{176 176 176}
{224 224 224}
{33 33 33}
{81 81 81}
{129 129 129}
{177 177 177}
{225 225 225}
{34 34 34}
{82 82 82}
{130 130 130}
{178 178 178}
{226 226 226}
{35 35 35}
{83 83 83}
{131 131 131}
{179 179 179}
{227 227 227}
{36 36 36}
{84 84 84}
{132 132 132}
{180 180 180}
{228 228 228}
{37 37 37}
{85 85 85}
{133 133 133}
{181 181 181}
{229 229 229}
{38 38 38}
{86 86 86}
{134 134 134}
{182 182 182}
{230 230 230}
{39 39 39}
{87 87 87}
{135 135 135}
{183 183 183}
{231 231 231}
{40 40 40}
{88 88 88}
{136 136 136}
{184 184 184}
{232 232 232}
{41 41 41}
{89 89 89}
{137 137 137}
{185 185 185}
{233 233 233}
{42 42 42}
{90 90 90}
{138 138 138}
{186 186 186}
{234 234 234}
{43 43 43}
{91 91 91}
{139 139 139}
{187 187 187}
{235 235 235}
{44 44 44}
{92 92 92}
{140 140 140}
{188 188 188}
{236 236 236}
{45 45 45}
{93 93 93}
{141 141 141}
{189 189 189}
{237 237 237}
{46 46 46}
{94 94 94}
{142 142 142}
{190 190 190}
{238 238 238}
{47 47 47}
{95 95 95}
{143 143 143}
{191 191 191}
{239 239 239}}}
{ "Gray2" { { 50 50 50 } { 81 81 81 } { 112 112 112 } { 143 143 143 } { 174 174 174 } { 205 205 205 } } }
{ "Gray3" { { 80 80 80 } { 99 99 99 } { 118 118 118 } { 137 137 137 } { 156 156 156 } { 175 175 175 } } }

	}
    }
    
    method getMaps {} {
	set maps {}
	for {set i 0} { $i < [llength $colorMaps]} {incr i} {
	    lappend maps [list [lindex [lindex $colorMaps $i] 0] $i]
	}
	return [list $maps]
    }
    
    method ui {} { 
      global $this-resolution
      
      set w .ui[modname]
      
      if {[winfo exists $w]} { 
          return
      } 
      
      set type ""
      
      sci_toplevel $w 
      wm minsize $w 200 50 
      
      sci_frame $w.f -relief flat -borderwidth 2
      pack $w.f -side top -expand yes -fill x 
      
      sci_frame $w.f.f1 -relief sunken -height 40  -borderwidth 2 
      pack $w.f.f1 -side right -padx 2 -pady 2 -expand yes -fill x


      canvas $w.f.f1.canvas -bg "#ffffff" -height 40 
      pack $w.f.f1.canvas -anchor w -expand yes -fill x

      TooltipMultiline $w.f.f1.canvas \
          "The red line represents the alpha value.  Use the left mouse button to add a\n" \
          "node for editing the line or to move an existing node.  You can use the\n" \
          "right mouse button to delete a node.  Alpha defaults to 0.5."

      sci_label $w.l0 -text "Click above to adjust alpha."
      pack $w.l0 -anchor c
      
      sci_frame $w.f3 -relief groove -borderwidth 2
      pack $w.f3 -side top -anchor c -expand yes -fill x -padx 2
      sci_scale $w.f3.s -orient horizontal -from -1 -to 1 -showvalue true \
          -label "Shift" -variable $this-gamma -resolution 0.01 -tickinterval 1
      pack $w.f3.s -expand yes -fill x -padx 2
      
      Tooltip $w.f3.s "Skews the color map to the left or right."

      sci_scale $w.f3.s2 -from 2 -to 256 -state normal \
        -orient horizontal  -variable $this-resolution -label "Resolution"
      pack $w.f3.s2 -expand yes -fill x -pady 2 -padx 2

      Tooltip $w.f3.s2 "Sets the number of unique colors used in the color map."
      
      bind $w.f3.s2 <ButtonRelease> "$this update; $this-c needexecute"

      sci_frame $w.f2 -relief groove -borderwidth 2
      pack $w.f2 -padx 2 -pady 2 -expand yes -fill both
      
      make_labeled_radio $w.f2.types "ColorMaps" "$this change" \
          top 2 $this-mapName [getColorMapNames]
          
      pack $w.f2.types -expand yes -fill both

      sci_frame $w.f4 -relief groove -borderwidth 2
      pack $w.f4 -padx 2 -pady 2 -expand yes -fill x

      sci_checkbutton $w.f4.faux -text "Opacity Modulation (Faux Shading)" -relief flat \
                -variable $this-faux -onvalue 1 -offvalue 0 \
                -anchor w -command "$this-c needexecute"
            pack $w.f4.faux -side top -fill x -padx 4
      Tooltip $w.f4.faux "Modulates color components based on the given opacity curve."

      sci_checkbutton $w.f4.reverse -text "Reverse the colormap" -relief flat \
                -variable $this-reverse -onvalue 1 -offvalue 0 \
                -anchor w -command "$this change"
            pack $w.f4.reverse -side top -fill x -padx 4
      Tooltip $w.f4.reverse "Reverse the colormap (not the alpha)"

      bind $w.f.f1.canvas <Expose> "$this canvasExpose"
      bind $w.f.f1.canvas <Button-1> "$this selectNode %x %y"
      bind $w.f.f1.canvas <B1-Motion> "$this moveNode %x %y"
      bind $w.f.f1.canvas <Button-3> "$this deleteNode %x %y"
      bind $w.f.f1.canvas <ButtonRelease> "$this update; $this-c needexecute"
      bind $w.f3.s <ButtonRelease> "$this change"
      $this update
      
      set cw [winfo width $w.f.f1.canvas]

      makeSciButtonPanel $w $w $this
      moveToCursor $w
    }
   
    method change {} {
	$this update
	$this-c needexecute
    }

    method selectNode { x y } {
	set w .ui[modname]
	set c $w.f.f1.canvas
	set curX $x
	set curY $y
	
	set selected [$c find withtag current]
	set index [lsearch [set $this-nodeList] $selected]
	if { $index == -1 } {
	    makeNode $x $y
	    set index [lsearch [set $this-nodeList] $selected]
	} 
	set loc [$c coords $selected]
	if { $loc != "" } {
	    set curX [expr ([lindex $loc 0]+[lindex $loc 2])*0.5]
	    set curY [expr ([lindex $loc 1]+[lindex $loc 3])*0.5]
	}
    }
	

    method makeNode { x y } {
	set w .ui[modname]
	set c $w.f.f1.canvas
	set new [$c create oval [expr $x - 5] [expr $y - 5] \
		     [expr $x+5] [expr $y+5] -outline white \
		     -fill black -tags node]
	set selected $new
	$this nodeInsert $new [list $x $y] $x
	$this drawLines
    }

    method drawLines { } {
	global $this-nodeList
	global $this-positionList
	set w .ui[modname]
	set canvas $w.f.f1.canvas
	$canvas delete line
	set cw [winfo width $canvas]
	set ch [winfo height $canvas]
	set x 0
	set y [expr $ch/2]
	for { set i 0 } { $i < [llength [set $this-nodeList]]} {incr i} {
	    set p [lindex [set $this-positionList] $i]
	    $canvas create line $x $y [lindex $p 0] [lindex $p 1] \
		-tags line -fill red
	    set x [lindex $p 0]
	    set y [lindex $p 1]
	}
	$canvas create line $x $y $cw [expr $ch/2]  -tags line -fill red
	$canvas raise node line
    }

    method drawNodes { } {
	global $this-nodeList
	global $this-positionList
	set w .ui[modname]
	set c $w.f.f1.canvas
	set cw [winfo width $c]
	set ch [winfo height $c]
	for {set i 0} { $i < [llength [set $this-nodeList]] } {incr i} {
	    set x [lindex [lindex [set $this-positionList] $i] 0 ]
	    set y [lindex [lindex [set $this-positionList] $i] 1 ]
	    if { $x < 0 } {
		set x 0
	    } elseif { $x > $cw } {
		set x $cw 
	    }

	    if { $y < 0 } {
		set y 0
	    } elseif { $y > $ch } {
		set y $ch 
	    }
	    set new [$c create oval [expr $x - 5] [expr $y - 5] \
		     [expr $x+5] [expr $y+5] -outline white \
		     -fill black -tags node]
	    set $this-nodeList [lreplace [set $this-nodeList] $i $i $new]
	}
    }

    method nodeInsert { n p x} {
	global $this-nodeList
	global $this-positionList
	set index 0
	for { set i 0 } { $i < [llength [set $this-nodeList]] } { incr i } { 
	    if { $x < [lindex [lindex [set $this-positionList] $i] 0 ] } {
		break;
	    } else {
		incr index
	    }
	}
	set $this-nodeList [linsert [set $this-nodeList] $index $n]
	set $this-positionList [linsert [set $this-positionList] $index $p]
    }
    
    method moveNode { x y } {
	global $this-nodeList
	global $this-positionList
	set w .ui[modname]
	set c $w.f.f1.canvas
	set cw [winfo width $c]
	set ch [winfo height $c]
	if { $curX + $x-$curX < 0 } { set x 0 }
	if { $curX + $x-$curX > $cw } { set x $cw}
	if { $curY + $y-$curY < 0 } { set y 0 }
	if { $curY + $y-$curY > $ch } { set y $ch }
	set i [lsearch  [set $this-nodeList] $selected]
	if { $i != -1 } {
	    set l [lindex [set $this-nodeList] $i]
	    $c move $l [expr $x-$curX] [expr $y-$curY]
	    set curX $x
	    set curY $y
	    set $this-nodeList [lreplace [set $this-nodeList] $i $i]
	    set $this-positionList [lreplace [set $this-positionList] $i $i]
	    nodeInsert $l [list $x $y] $x
	    $this drawLines
	}
    }

    method lreverse { stuff } {
	set size [llength $stuff]
	set result {}
	for {set i 0}  {$i < $size} {incr i} {
	    set result [concat [list [lindex $stuff $i]] $result]
	}
	return $result
    }

    method findByName_aux { cname } {
	set size [llength $colorMaps]
	for {set i 0}  {$i < $size} {incr i} {
	    set cname1 [lindex [lindex $colorMaps $i] 0]
	    if {$cname == $cname1} { return $i }
	}
	return 0
    }

    method findByName { cname } {
	set index [findByName_aux $cname]
	set color [lindex [lindex $colorMaps $index] 1]
	if {[set $this-reverse]} { set color [lreverse $color] }
	return $color
    }

    method getColorMapNames {} {
	set size [llength $colorMaps]
	set result {}
	for {set i 0}  {$i < $size} {incr i} {
	    set cname [lindex [lindex $colorMaps $i] 0]
	    set result [concat $result [list [list $cname $cname]]]
	}
	return $result
    }
	
    method lookupOldIndex {a b c} {
	global $this-mapType
	set index [set $this-mapType]
	# Old name, new name, reverse?  Note that the order these are
	# listed is important because mapType is an index into this
	# list.
	set remap {
	    { "Gray" "Gray" 0 }
	    { "Inverse Gray" "Gray" 1 }
	    { "Rainbow" "Old Rainbow" 1 }
	    { "Inverse Rainbow " "Old Rainbow" 0 }
	    { "Darkhue" "Darkhue" 0 }
	    { "Inverse Darkhue" "Darkhue" 1 }
	    { "Lighthue" "Lighthue" 0 }
	    { "Blackbody" "Blackbody" 0 }
	    { "Inverse Blackbody" "Blackbody" 1 }
	    { "Don" "Don" 0 }
	    { "Dark Gray" "Dark Gray" 0 }
	    { "Red Tint" "Red Tint" 0 }
	    { "Orange Tint" "Orange Tint" 0 }
	    { "Yellow Tint" "Yellow Tint" 0 }
	    { "Green Tint" "Green Tint" 0 }
	    { "Blue Tint" "Blue Tint" 0 }
	    { "Purple Tint" "Purple Tint" 0 }
	    { "BP Seismic" "BP Seismic" 0}
           }
	set entry [lindex $remap $index]
	set $this-mapName [lindex $entry 1]
	set $this-reverse [lindex $entry 2]
    }

    method deleteNode { x y } {
	global $this-nodeList
	global $this-positionList
	set w .ui[modname]
	set c $w.f.f1.canvas
	set l [$c find withtag current]
	set i [lsearch  [set $this-nodeList] $l]
	if { $i != -1 } {
	    $c delete current
	    set $this-nodeList [lreplace [set $this-nodeList] $i $i]
	    set $this-positionList [lreplace [set $this-positionList] $i $i]
	    $this drawLines
	}
    }
    
    method update { } {
      $this SetColorMap
      $this redraw
      set selected -1
    }
    
    method getColorMapString { } {
      set colors [findByName [set $this-mapName]]
      set csize [llength $colors]
      set scolors [join $colors]

      global $this-positionList
      global $this-width
      global $this-height
      set cw [set $this-width]
      set ch [set $this-height]
      set alphas [join [set $this-positionList]]
      set asize [llength [set $this-positionList]]

      return "$csize $scolors $asize $cw $ch $alphas"
    }

    method close {} {
      set w .ui[modname]
      set exposed 0
      destroy $w
    }
    
    method canvasExpose {} {
	set w .ui[modname]
	
	if { [winfo viewable $w.f.f1.canvas] } { 
	    if { $exposed } {
		return
	    } else {
		set exposed 1
		$this drawNodes
		$this drawLines
		$this redraw
	    } 
	} else {
	    return
	}
    }
    
    method redraw {} {
      global $this-width
      global $this-height
      set w .ui[modname]
      
      if {![winfo exists $w]} { 
          return
      } 

      set n [llength $colorMap]
      set canvas $w.f.f1.canvas
      $canvas delete map
      set cw [winfo width $canvas]
      set $this-width $cw
      set ch [winfo height $canvas]
      set $this-height $ch
      set dx [expr $cw/double($n)] 
      set x 0
      for {set i 0} {$i < $n} {incr i 1} {
          set color [lindex $colorMap $i]
          set r [lindex $color 0]
          set g [lindex $color 1]
          set b [lindex $color 2]
          set c [format "#%02x%02x%02x" $r $g $b]
          set oldx $x
          set x [expr ($i+1)*$dx]
          $canvas create rectangle \
            $oldx 0 $x $ch -fill $c -outline $c -tags map
      }
      set taglist [$canvas gettags all]
      set i [lsearch $taglist line]
      if { $i != -1 } {
          $canvas lower map line
      }
    }

    method SetColorMap {} {
      global $this-resolution
      global $this-mapName
      set colorMap {}
      set map [findByName [set $this-mapName]]
      set currentMap {}
      set currentMap [$this makeNewMap $map]
      set n [llength $currentMap]

      set res [set $this-resolution]

      set frac [expr ($n-1)/double($res-1)]
      for { set i 0 } { $i < $res  } { incr i} {
        if { $i == 0 } {
          set color [lindex $currentMap 0]
          lappend color [$this getAlpha $i]
        } elseif { $i == [expr ($res - 1)] } {
          set color [lindex $currentMap [expr ($n - 1)]]
          lappend color [$this getAlpha $i]
        } else {
          set index [expr int($i * $frac)]
          set t [expr ($i * $frac)-$index]
          set c1 [lindex $currentMap $index]
          set c2 [lindex $currentMap [expr $index + 1]]
          set color {}
          for { set j 0} { $j < 3 } { incr j} {
            set v1 [lindex $c1 $j]
            set v2 [lindex $c2 $j]
            lappend color [expr int($v1 + $t*($v2 - $v1))]
          }
          lappend color [$this getAlpha $i]
        }
        lappend colorMap $color
      }
    }
    

    method makeNewMap { currentMap } {
	global $this-gamma
	global $this-r

	set res [set $this-resolution]
	set newMap {}
	set m [expr int($res + abs( [set $this-gamma] )*(255 - $res))]
	set n [llength $currentMap]
	if { $m < $n } { set m $n }
	set frac [expr double($n-1)/double($m - 1)]
	for { set i 0 } { $i < $m  } { incr i} {
	    if { $i == 0 } {
		set color [lindex $currentMap 0]
	    } elseif { $i == [expr ($m -1)] } {
		set color [lindex $currentMap [expr ($n - 1)]]
	    } else {
		set index_double [$this modify [expr $i * $frac] [expr $n-1]]
		
		set index [expr int($index_double)]
		set t  [expr $index_double - $index]
		set c1 [lindex $currentMap $index]
		set c2 [lindex $currentMap [expr $index + 1]]
		set color {}
		for { set j 0} { $j < 3 } { incr j} {
		    set v1 [lindex $c1 $j]
		    set v2 [lindex $c2 $j]
		    lappend color [expr int($v1 + $t*($v2 - $v1))]
		}
	    }
	    lappend newMap $color
	}
	return $newMap
    }
    
    method modify {  i range } {
      global $this-gamma
      
      set val [expr $i/double($range)]
      set bp [expr tan( 1.570796327*(0.5 + [set $this-gamma]*0.49999))]
      set index [expr pow($val,$bp)]
      return $index*$range
    }
    
    method getAlpha { index } {
      global nodeList
      global $this-positionList
      global $this-width
      global $this-height

      set cw [set $this-width]
      set ch [set $this-height]
      set m [llength [set $this-nodeList]]
      set dx [expr $cw/double([set $this-resolution])] 
      set xpos [expr int($index*$dx) + 0.5 * $dx]
      set x 0.0
      set y $ch/2.0
      set i 0
      for {set i 0} {$i <= $m} {incr i 1} {
        set newx 0
        set newy $ch/2.0
        if { $i != $m } {
          set newx [lindex [lindex [set $this-positionList] $i] 0]
          set newy [lindex [lindex [set $this-positionList] $i] 1]
        } else {
          set newy $ch/2.0
          set newx $cw
        }

        if { $xpos < $newx } {
          set frac [expr ( $xpos - $x )/double( $newx - $x )]
          return [expr 1.0 -  (($y + $frac * ($newy - $y))/double($ch))]
        }
        set x $newx
        set y $newy
      }
    }
}

