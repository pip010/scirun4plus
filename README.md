SCIRun4.7+
===========

SCIRun (4.7) plus additional modules

[![DOI](https://zenodo.org/badge/15745408.svg)](https://zenodo.org/badge/latestdoi/15745408)

Modules
=======

SolveBiotSavartContour : Given a domain mesh and coil mesh in the for of curve-mesh the module calculates the 
B-field and A-field per mesh node using BiotSvart method for close contours. 
Piece-wise integration is conducted per domain mesh node and the constant value per coil model edge 
is considered as the current through the coil at any certain moment. Positive current is assumed to designate 
integration in direction as the order of edge nodes, while a negative one in the opposite direction.
File(s): 
/src/Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.cc/.h
/src/DataFlow/XML/SolveBiotSavartContour.xml
/src/DataFlow/Modules/Math/SolveBiotSavartContour.cc

ModelCoil : The module serves as a helper for generating coil mesh in accordance of predefined coil(s) configuration.
File(s): 
/src/Core/Algorithms/Math/TMS/ModelGenericCoilAlgo.cc/.h
/src/DataFlow/GUI/ModelTMSCoilSingle.tcl
/src/DataFlow/XML/ModelTMSCoilSingle.xml
/src/DataFlow/Modules/TMS/ModelTMSCoilSingle.cc
/src/DataFlow/GUI/ModelTMSCoilSpiral.tcl
/src/DataFlow/XML/ModelTMSCoilSpiral.xml
/src/DataFlow/Modules/TMS/ModelTMSCoilSpiral.cc
/src/DataFlow/GUI/ModelTMSCoilDipole.tcl
/src/DataFlow/XML/ModelTMSCoilDipole.xml
/src/DataFlow/Modules/TMS/ModelTMSCoilDipole.cc


RunExtProcess : This module allows executing external process and perform basic I/O via system files or standart I/O stream.
File:NA

Utility
=====

CoilModelCalibration: Utility app that reports the magnetic field magnitude along the central axis through circular coil. 
Reported values are at 1,2,3,5,8 cm away from the circle center along the central axis. 
 
 Usage: ./CoilModelCalibration type lod istep radius
	 Reports magnetic field magnitude at axial distance of 1, 2, 3, 5, 8 cm.
	 where, type: 1 - circular coil  
	              2 - spiral coil 
	              3 - dipoles model 
	              0 - analytical 
	 lod: level of detail for coil geom [integer] 
	 istep: integration step [real]
	 radius: radius of the coil [real]

Files(s):
/src/Applications/Utils/CoilModelCalibration.cc

Build
=====
 ./build.sh (builds in release using only 1 core on your CPU)
 ./build.sh -h (for help)
 ./build.sh --debug -j4 (builds in debug on 4 cores CPU)
 ./build.sh --hybrid (builds in release with pdb debug symbols, handy for profiling and debug release only issues)
 
*PreReq. Linux : http://scirundocwiki.sci.utah.edu/SCIRunDocs/index.php/Ubuntu_Prerequisites_%28Packages_to_Install%29

Ubuntu
=====

Ubuntu comes with a relatively minimal install, so several packages are required.
Ubuntu 12.x packages

Tested on Ubuntu 12.04. It is easy to get the necessary packages listed below via the Synaptic package manager or the command line instructions below:

sudo apt-get install subversion cmake-qt-gui cmake-curses-gui build-essential \
                     libxft-dev libxi-dev libxmu-headers freeglut3-dev \
                     libtiff4-dev

Ubuntu 11.x packages

Tested on Ubuntu 11.10 and 11.04. It is easy to get the necessary packages listed below via the Synaptic package manager or the command line instructions below:

sudo apt-get install subversion cmake-qt-gui cmake-curses-gui build-essential \
                     libxft-dev libxi-dev libxmu-headers freeglut3-dev \
                     libtiff4-dev

Ubuntu 10.x packages

Tested on Ubuntu 10.04 and 10.10. It is easy to get the necessary packages listed below via the Synaptic package manager or the command line instructions below:

sudo apt-get install subversion cmake-qt-gui cmake-curses-gui build-essential \
                     libxft-dev libxi-dev libxmu-headers libglut3-dev \
                     libtiff4-dev

Tips
====

Although, SCIRun 4.x is multi-OS in principle. I have experience only building it on Linux. 
Besides the dependencies you need to resolve first pay attention to use CMake 2.8 or higher and
to build in the home folder. Don't try to install !
