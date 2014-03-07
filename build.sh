#!/bin/bash
#  
#  For more information, please see: http://software.sci.utah.edu
#  
#  The MIT License
#  
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
#  University of Utah.
#  
#  License for the specific language governing rights and limitations under
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
#    File   : build.sh
#    Author : McKay Davis
#    Date   : Tue Dec 19 15:35:50 2006


# SCIRun build script
#
# This script will build SCIRun from scratch
#

CMAKE_MAJOR_VERSION=2
CMAKE_MINOR_VERSION=8
CMAKE_PATCH_VERSION=11
CMAKE_VERSION="${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}"
CMAKE="cmake-$CMAKE_VERSION"

OSX_TARGET_VERSION="10.5"
TARGET_VERSION_SDK="/Developer/SDKs/MacOSX10.5.sdk"

printhelp() {
    echo -e "--cmake-args=[string]\tCommand line arguments to CMake"
    echo -e "--debug\t\t\tBuilds SCIRun with debug symbols"
    echo -e "--release\t\tBuilds SCIRun without debug symbols (default)"
    echo -e "--without-itk\t\tTurns off building ITK and the Insight SCIRun package"
    echo -e "--without-meshing-pipeline\tDownloads necessary thirdparties and builds meshing pipeline"
    echo -e "--without-tetgen\tTurns off tetgen support. By default, tetgen is built with SCIRun"
    echo -e "--without-lapack\t\tTurn on SCIRun LAPACK support."
    echo -e "--get-cmake\t\tDownloads and builds $CMAKE, regarless of the version installed on the system"
    echo -e "--static\t\tBuilds static libs"
    echo -e "--packages=[string]\tComma seperated string specifies which SCIRun packages are built"
    echo -e "--no-ctest\t\tDo not run ctest"
    echo -e "--64bit\t\t\tForce a 64bit SCIRun version [OSX only]"
    echo -e "--set-osx-version-min\tTarget a minimum Mac OSX version (currently $OSX_TARGET_VERSION) [OSX only]"
    echo -e "--universal\t\tEnable universal build for Mac OSX (NB: OSX 10.6 and later do not support ppc)"
    echo -e "--documentation\t\tEnable building documentation (requires LaTeX)"
    echo -e "-j#\t\t\tRuns # parallel make processes when building"
    echo -e "-?\t\t\tThis help"
    exit 0
}

# will cause the script to bailout if the passed in command fails
try () {
  $*
  if [[ $? != "0" ]]; then
      echo -e "\n***ERROR in build script\nThe failed command was:\n$*\n"
      exit 1
  fi
}

trybuild () {
  $*
  if [[ $? != "0" ]]; then
      echo -e "Building SCIRun returned an error\n"
      echo -e "Either the code failed to build properly or\n"
      echo -e "the testing programs failed to complete without\n"
      echo -e "every single test program passing the test.\n"
      exit 1
  fi
}
# functionally equivalent to try(),
# but it prints a different error message
ensure () {
  $* >& /dev/null
  if [[ $? != "0" ]]; then
      echo -e "\n***ERROR, $* is required but not found on this system\n"
      exit 1
  fi
}

# Try to find a version of cmake
find_cmake() {
    if [[ "$getcmake" != "1" ]]; then
        cmakebin=`which cmake`
        ctestbin=`which ctest`
    fi

    download=0
    #if it is not found
    if [[ ! -e "$cmakebin" ]]; then
        download=1
    else
    # see if cmake is up-to-date
        version=`$cmakebin --version | cut -d ' ' -f 3 | sed -e "s/[^[:digit:].]//g"`
        echo "$cmakebin version $version found"
        major_version=${version:0:1}
        minor_version=${version:2:1}
        if [[ $major_version -le $CMAKE_MAJOR_VERSION && $minor_version -lt $CMAKE_MINOR_VERSION ]] ; then
            download=1
        fi
    fi

    if [[ $download -eq 1 ]]; then
    # then look for our own copy made by this script previously
        cmakebin=$DIR/cmake/local/bin/cmake
        ctestbin=$DIR/cmake/local/bin/ctest
        try mkdir -p $DIR/cmake/
        try cd $DIR/cmake

        if [[ -e "$cmakebin" ]]; then
        # see if local cmake install is compatible
            version=`$cmakebin --version | cut -d ' ' -f 3 | sed -e "s/[^[:digit:].]//g"`
            echo "$cmakebin version $version found"
            major_version=${version:0:1}
            minor_version=${version:2:1}
            if [[ $major_version -ge ${CMAKE_MAJOR_VERSION} && $minor_version -ge ${CMAKE_MINOR_VERSION} ]] ; then
                download=0
            fi
        fi

        if [[ $download -eq 1 ]]; then
        # try to download and build our own copy in local
            try $getcommand http://www.cmake.org/files/v${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}/${CMAKE}.tar.gz
            try tar xvzf ${CMAKE}.tar.gz
            try cd $CMAKE
            try ./bootstrap --prefix=$DIR/cmake/local
            try make $makeflags
            try make install
        fi
    fi

    ctestbin="$ctestbin --ctest-config=${DIR}/src/CTestConfig.cmake"

    echo cmakebin=$cmakebin
    echo ctestbin=$ctestbin
    ensure $cmakebin --version
}

find_svn() {
    # Try to find a version of svn
    svnbin=`which svn`
    echo svnbin=$svnbin
}

update_scirun_src() {
    echo "Updating SCIRUN source..."

    # If the source tree is from SVN, update it too
    #if [[ -e "$DIR/src/.svn" ]] && [[ -e "$svnbin" ]]; then 
    #    try svn update "$DIR/src"
    #fi
}
    
configure_scirun_bin() {
    try cd $DIR/bin
    
    if [[ "$builddataflow" == "1" ]]; then
      buildsharedlibs="1"
    fi

    # Put common SCIRun build options here:
    local COMMON_BUILD_OPTS="-DBUILD_SHARED_LIBS=$buildsharedlibs \
                             -DBUILD_DATAFLOW=$builddataflow \
                             -DCMAKE_BUILD_TYPE=$buildtype \
                             -DLOAD_PACKAGE=$scirunpackages \
                             -DWITH_X11=$withx11 \
                             -DBUILD_TESTING=$testnetworks \
                             -DBUILD_UTILS=$buildutils \
                             -DBUILD_BIOMESH3D=$buildmeshing \
                             -DWITH_TETGEN=$buildtetgen \
                             -DWITH_LAPACK=$buildlapack \
                             -DWITH_ITK=$builditk \
                             -DCMAKE_VERBOSE_MAKEFILE=$verbosebuild"

    if [[ "$osx" == "1" ]]; then
      if [[ "$setosxmin" == "1" ]]; then
        # set up target SDK for (useful for building binary installers)
        COMMON_BUILD_OPTS="$COMMON_BUILD_OPTS -DCMAKE_OSX_DEPLOYMENT_TARGET=$OSX_TARGET_VERSION -DCMAKE_OSX_SYSROOT=$TARGET_VERSION_SDK"
      fi

      # mac architecture-specific cmake options
      if [[ "$universal" == "1" ]]; then
        # explicitly set architecture types for universal builds
        if [[ "$darwin64" == "1" ]]; then
          $cmakebin ../src -DSCIRUN_ENABLE_64BIT=ON -DCMAKE_OSX_ARCHITECTURES="ppc64;x86_64" -DCMAKE_C_FLAGS="-arch x86_64 -arch ppc64 " -DCMAKE_CXX_FLAGS="-arch x86_64 -arch ppc64 " -DCMAKE_EXE_LINKER_FLAGS="-arch x86_64 -arch ppc64 " $COMMON_BUILD_OPTS $cmakeargs
        else
          $cmakebin ../src -DCMAKE_OSX_ARCHITECTURES="ppc;i386" -DCMAKE_C_FLAGS="-arch i386 -arch ppc " -DCMAKE_CXX_FLAGS="-arch i386 -arch ppc " -DCMAKE_EXE_LINKER_FLAGS:="-arch i386 -arch ppc " $COMMON_BUILD_OPTS $cmakeargs
        fi
      else
        if [[ "$darwin64" == "1" ]]; then
          $cmakebin ../src -DSCIRUN_ENABLE_64BIT=ON -DCMAKE_OSX_ARCHITECTURES="x86_64" -DCMAKE_C_FLAGS="-m64 " -DCMAKE_CXX_FLAGS="-m64 " -DCMAKE_EXE_LINKER_FLAGS="-m64 " $COMMON_BUILD_OPTS $cmakeargs
        elif [[ "$darwin32" == "1" ]]; then
          # explicitly force a 32bit build
          $cmakebin ../src -DCMAKE_OSX_ARCHITECTURES="i386" -DCMAKE_C_FLAGS="-m32 " -DCMAKE_CXX_FLAGS="-m32 " -DCMAKE_EXE_LINKER_FLAGS:="-m32 " $COMMON_BUILD_OPTS $cmakeargs
        else
          # use default system settings
          $cmakebin ../src $COMMON_BUILD_OPTS $cmakeargs
        fi
      fi

    else
      # all other *nix builds
      $cmakebin ../src $COMMON_BUILD_OPTS $cmakeargs
    fi
}

build_scirun_bin() {
    if [[ "$usectest" == "1" ]] && [[ -e "$ctestbin" ]]; then
        echo "Building SCIRun using ctest..."
        trybuild $ctestbin -VV -D Experimental -A $DIR/bin/CMakeCache.txt
    else 
	echo "Building SCIRun using make..."
        trybuild make $makeflags
    fi
}


######### build.sh script execution starts here

export DIR=`pwd`
export SCIRUN_THIRDPARTY_DIR=$DIR/bin

linux=0
osx=0

if test `uname` = "Darwin"; then
    getcommand="curl -OL"
    osx=1
elif test `uname` = "Linux"; then
    getcommand="wget"
    linux=1
else
    echo "Unsupported system.  Please run on OSX or Linux"
    exit 1
fi

buildtype="Release"
buildsharedlibs="1"
builddataflow="1"
makeflags=""
cmakeflags=""
builditk="1"
buildmeshing="1"
buildtetgen="1"
buildlapack="1"
getcmake=0
withx11="1"
scirunpackages="SCIRun,BioPSE,Teem,MatlabInterface"
cmakeargs=""
usectest="0"
testnetworks="0"
buildutils="1"
darwin64="0"
darwin32="0"
setosxmin="0"
verbosebuild="0"
universal="0"
documentation="0"

echo "Parsing arguments..."
while [[ "$1" != "" ]]; do
    case "$1" in
        --universal)
            universal="1";;
        --documentation)
            documentation="1";;
        --debug)
            buildtype="Debug";;
        --release)
            buildtype="Release";;
        --without-itk) 
            builditk="0";;
        --without-tetgen)
            buildtetgen="0";;
        --without-lapack)
            buildlapack="0";;
        --without-meshing-pipeline)
            buildmeshing="0";;
        --get-cmake) 
            getcmake="1";;
        --no-ctest) 
            usectest="0";;
        --ctest) 
            usectest="1";;
        --test-networks)
            testnetworks="1";;
        --static) 
            buildsharedlibs="0";;
       --64bit)
            if [[ "$osx" == "1" ]]; then
              darwin64="1"
            else
              echo "Only OSX supporsts the selection between 32 and 64bit architectures."
              echo "On other platforms SCIRun uses the architecture of the operating system."
            fi;;
       --32bit)
            if [[ "$osx" == "1" ]]; then
              darwin32="1"
            else
              echo "Only OSX supporsts the selection between 32 and 64bit architectures."
              echo "On other platforms SCIRun uses the architecture of the operating system."
            fi;;
       --set-osx-version-min)
            if [[ "$osx" == "1" ]]; then
              setosxmin="1"
            else
              echo "Only OSX supports the --set-osx-version-min flag."
            fi;;
       --verbose)
            verbosebuild="1";;
       --packages=*)
            scirunpackages=`echo $1 | awk -F '\=' '{print $2}'`;;
        --cmake-args=*)
            cmakeargs=`echo $1 | cut -c 14-`;;
        -j*) 
            makeflags="${makeflags} $1";;
        -D*) 
            cmakeflags="${cmakeflags} $1";;
        -?|--?|-help|--help) 
            printhelp;;
        *) 
            echo \`$1\' parameter ignored;;
    esac
    shift 1
done

cmakeargs="${cmakeargs} ${cmakeflags}"

if [[ "$documentation" != "0" ]] ; then
  cmakeargs="${cmakeargs} -DBUILD_DOCUMENTATION=ON"
fi

echo "CMake args: $cmakeargs"
echo "Build Type: $buildtype"
echo "Build SCIRun Dataflow: $builddataflow"
if [[ "$builddataflow" != "0" ]]; then
    echo "SCIRun Packages: $scirunpackages"
fi
echo "Build Insight ToolKit: $builditk"
echo "Get CMake: $getcmake"
echo "Make Flags: $makeflags"
echo "Build Shared Libs: $buildsharedlibs"
echo "With X11: $withx11"
echo "With LAPACK: $buildlapack"

ensure make --version # ensure make is on the system

find_cmake

find_svn

update_scirun_src

configure_scirun_bin

build_scirun_bin
