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



SCIRUN_VERSION=3.2.0
CMAKE=cmake-2.4.5

printhelp() {
    echo -e "--nodashboard\t\tDo not commit results to dashboard"
    echo -e "--newimages\t\tRender a new set of baseline images"
    echo -e "--imagedir\tSet the directory with the baseline images"
    echo -e "--datadir\tSet the directory with scirun data"
    echo -e "-?\t\t\tThis help"
    exit 0
}

# will cause the script to bailout if the passed in command fails
try () {
  $*
  if [ $? != "0" ]; then
      echo -e "\n***ERROR in build script\nThe failed command was:\n$*\n"
      exit 1
  fi
}

trybuild () {
  $*
  if [ $? != "0" ]; then
      echo -e "Building SCIRun returned an error\n"
      echo -e "Either the code failed to build properly or\n"
      echo -e "the testing programs failed to complete without\n"
      echo -e "every single test program passing the test.\n"
  fi
}
# functionally equivalent to try(),
# but it prints a different error message
ensure () {
  $* >& /dev/null
  if [ $? != "0" ]; then
      echo -e "\n***ERROR, $* is required but not found on this system\n"
      exit 1
  fi
}




# Try to find a version of cmake
find_cmake() {
  cmakebin=`which cmake`
  ctestbin=`which ctest`
}

do_ctest() {


  if [ "$newimages" == "1" ]; then
    rm -f $DIR/Testing/*.png
    trybuild $cmakebin . $sitename -DWITH_WXWIDGETS=OFF -DBUILD_TESTING=ON -DRUN_UNIT_TESTS=ON -DRUN_CLASS_TESTS=ON -DRUN_SAMPLE_TESTS=ON -DTEST_IMAGES_DIR=$DIR/Testing -DCOMPARE_IMAGES=OFF -DBUILD_SEG3D=OFF
#    trybuild make -j4
    trybuild $ctestbin $submit
    rm -rf $imagedir
    mkdir $imagedir
    cp $DIR/Testing/*.png $imagedir
  else
    trybuild $cmakebin . $sitename -DWITH_WXWIDGETS=OFF -DBUILD_TESTING=ON -DRUN_UNIT_TESTS=ON -DRUN_CLASS_TESTS=ON -DRUN_SAMPLE_TESTS=ON -DTEST_IMAGES_DIR=$DIR/Testing -DCOMPARE_IMAGES=ON -DBASELINE_IMAGES_DIR=$imagedir -DBUILD_SEG3D=OFF
#    trybuild make -j4
    trybuild $ctestbin $submit  
  fi
}



######### build.sh script execution starts here

export DIR=`pwd`
linux=0
osx=0

if test `uname` = "Darwin"; then
    osx=1
elif test `uname` = "Linux"; then
    linux=1
else
    echo "Unsupported system.  Please run on OSX or Linux"
    exit 1
fi

submit="-D Experimental -ctest-config=$DIR/../src/CTestConfig.cmake -VV"
newimages=0

echo "Parsing arguments..."

imagedir=$DIR/Testing/BaselineImages/
export SCIRUN_DATA="$DIR/../../SCIRunData"

sitename="" 

while [ "$1" != "" ]; do
    case "$1" in
        --dashboard) 
          submit="-D Experimental -ctest-config=$DIR/../src/CTestConfig.cmake -VV";;
        --nightly)
          submit="-D Nightly -ctest-config=$DIR/../src/CTestConfig.cmake -VV";;
        --continuous)
          submit="-D Continuous -ctest-config=$DIR/../src/CTestConfig.cmake -VV";;          
        --nodashboard)
          submit="";;
        --newimages)
          newimages=1;;
        -?|--?|-help|--help) 
            printhelp;;
        --sitename=*)
          sitename=`echo $1 | awk -F '\=' '{print " -DSITE:STRING="$2}'`;;        
        --imagedir=*)
          imagedir=`echo $1 | awk -F '\=' '{print $2}'`;;
        --datadir=*)
          datadir=`echo $1 | awk -F '\=' '{print $2}'`
          export SCIRUN_DATA=`echo $1 | awk -F '\=' '{print $2}'`;;
        *) 
            echo \`$1\' parameter ignored;;
    esac
    shift 1
done

find_cmake

do_ctest 
