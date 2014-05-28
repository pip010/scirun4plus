scirun4plus
===========

scirun4 plus additional modules

modules
=======

SolveBiotSavartContour : Given a domain mesh and coil mesh in the for of curve-mesh the module calculates the B-field and A-field per mesh node using BiotSvart method for close contours. Piece-wise integration is conducted per domain mesh node and the constant value per coil model edge is considered as the current through the coil at any certain moment. Positive current is assumed to designate integration in direction as the order of edge nodes, while a negative one in the direction oposite the order of edge nodes.

ModelCoil : The module serves as a helper for generating coil mesh in accordance of predefined coil(s) configuration.

RunExtProcess : This module allows executing external process and perform basic I/O via system files or standart I/O stream.


build
=====
 ./build.sh (builds in release using only 1 core on your CPU)
 ./build.sh -h (for help)
 ./build.sh --debug -j4 (builds in debug on 4 cores CPU)
 ./build.sh --hybrid (builds in release with pdb debug symbols, handy for profiling and debug release only issues)
*PreReq. Linux : http://scirundocwiki.sci.utah.edu/SCIRunDocs/index.php/Ubuntu_Prerequisites_%28Packages_to_Install%29
