#!python
#
#  For more information, please see: http://software.sci.utah.edu
# 
# The MIT License
#
# Copyright (c) 2009 Scientific Computing and Imaging Institute,
# University of Utah.
#
#
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
#
#
# Tom Fogal, Tue Nov 11 12:50:30 MST 2008 
# Does a build of thirdparty and SCIRun trees.  Copies DLLs to the executable
# directory so that the tests run `out of the box'.
# Run with --help for usage information.

import logging
import os
from optparse import OptionParser
import shutil
import subprocess
import sys

#### Path Setup ####
# These should all be changed to be defaults modified by command line options,
# but I haven't made time for this yet.

# Root dir off of which everything is based.  Purely a convenience.
ROOT="C:/Users/admin/Documents"

# Where the branch is checked out to (or where it should be checked out to)
SCIRUN="%s/source/SCIRun-trunk" % (ROOT, )

# source directory within the checkout
SRC="%s/src" % (SCIRUN, )

# Path to the microsoft SDK that we should use for building.
SDK="C:/Program Files/Microsoft SDKs/Windows/v6.0A"

# TetGen binary location
TETGEN="%s/mesh_thirdparty/tetgen1.4.2/tetgen1.4.2/VS 2005/tetgen/x64/Release" % (ROOT)

#### End of path setup ####
# No machine-specific modification should be necessary below this line.

parser = OptionParser()
parser.add_option("-d", "--debug", dest="debug",
                  help="Enable debug mode.",
                  action="store_true", default=False)
parser.add_option("-t", "--no-thirdparty", dest="thirdparty",
                  help="Skip doing the thirdparty build.",
                  action="store_false", default=True)
parser.add_option("-c", "--no-cache", default=False,
                  help="Delete CMake cache before running cmake",
                  action="store_true", dest="del_cache")
parser.add_option("-C", "--no-cmake", default=True,
                  help="Skip running CMake",
                  action="store_false", dest="do_cmake")
parser.add_option("-r", "--run-tests", default=False,
                  help="Run the test suite after building.",
                  action="store_true", dest="run_tests")
parser.add_option("--clean", dest="clean",
                  help="Run a `make clean' step before building.",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

logging.basicConfig(level=logging.INFO)

build_type = "Release"
if options.debug:
    build_type = "Debug"

# The directory to build in.
# FIXME -- create the directory, or at least make sure it exists first!
BUILD="%s/build/trunk-%s" % (ROOT, build_type)

BUILDSTR="%s64" % (build_type, )

THIRDPARTY="%s/thirdparty.win64" % (SCIRUN)
BUILD_NAME="x86_64-windows-%s" % (build_type)
TEEM_SRC="%s/teem-1.9.0-src" % (THIRDPARTY)

bld_config = "%s|x64" % (build_type, )

print "Doing a", build_type, "in", BUILD
print "Build Configuration:", bld_config
print "With thirdparty from:", THIRDPARTY

def msbuild(solution, project):
    subprocess.check_call(["msbuild", "/m", "/nologo",
                           "/t:%s" % (project),
                           "/v:n", # changes the verbosity
                           "/maxcpucount:8",
                           "/p:Configuration=%s,Platform=x64" % (build_type),
                           solution], shell=True)

def build(solution, project):
    '''Builds via msbuild, but cleans first if the user wants.'''
    if options.clean:
        logging.info("Cleaning before build ...")
        msbuild(solution, "clean")
    logging.info("Building...")
    msbuild(solution, project)

def devenv(solution, project):
    subprocess.check_call(["devenv", solution,
                           "/build", bld_config,
                           "/project", project,
                           "/projectconfig", bld_config,
                           "/nologo"], bufsize=1, shell=True)
def rm_f(path):
    '''Removes a given path.  Non-fatal if it cannot be removed.'''
    try:
        os.unlink(path)
    except OSError:
        logging.warn("Could not delete '%s'" % (path))

# build tetgen via python?
# cl /c predicates.cxx
# cl /c tetgen.cxx
# link /out:tetgen.lib /machine:x64 predicates.obj tetgen.obj

def move_to_dir(path):
    if os.access(path, os.F_OK) is False:
        os.mkdir(path)
    os.chdir(path)
    print "Switched to:", os.getcwd()

if options.thirdparty:
    # Temporary fix -- don't build teem, because this code is for building teem
    # via CMake, and for today we're sticking to building teem through our
    # normal thirdparty VS `solution'.
#    move_to_dir("%s/build/teem-%s" % (ROOT, build_type))
#
#    if options.del_cache:
#        rm_f("CMakeCache.txt")
#
#        "-DBUILD_SHARED_LIBS:BOOL=OFF",
#    subprocess.check_call(["cmake", "-G", "Visual Studio 9 2008 Win64",
#        "-DBUILD_TYPE:STRING=%s" % (build_type),
#        "-DCMAKE_BUILD_TYPE:STRING=%s" % (build_type),
#        "-DCMAKE_VERBOSE_MAKEFILE:STRING=%s" % (build_type),
#        "-DCMAKE_INSTALL_PREFIX:PATH=%s/%s64" % (THIRDPARTY, build_type),
#        "%s" % (TEEM_SRC)])
#    build("TEEM.sln", "ALL_BUILD")
#    # Can't install from msbuild for some stupid reason.
#    devenv("TEEM.sln", "INSTALL")

    move_to_dir(THIRDPARTY)
    print "Switched to thirdparty dir:", os.getcwd()
    build("thirdparty.2008.sln", "finalize")

move_to_dir(BUILD)
print "Switched to build directory:", os.getcwd()
if options.del_cache:
    rm_f("CMakeCache.txt")

#        "-DBUILD_SHARED_LIBS:BOOL=ON",

if options.do_cmake:
    subprocess.check_call(["cmake", "-G", "Visual Studio 9 2008 Win64",
        '-DSCIRUN_THIRDPARTY_DIR:PATH=%s/%s' % (THIRDPARTY, BUILDSTR),
        "-DLOAD_PACKAGE:STRING=SCIRun,BioPSE,Teem,MatlabInterface",
        "-DCMAKE_CXX_FLAGS:STRING=-UHAVE_HASH_MAP /EHsc /MT",
        "-DSubversion_SVN_EXECUTABLE:PATH=", # purposely unset svn, for git
        "-DWITH_X11:BOOL=OFF",
        "-DBUILD_TYPE:STRING=%s" % (build_type),
        "-DCMAKE_BUILD_TYPE:STRING=%s" % (build_type),
        "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON",
        "-DBUILD_SEG3D:BOOL=OFF",
        "-DSCIRUN_64BITS_INDEX:BOOL=ON",
        "-DSCIRUN_ENABLE_64BIT:BOOL=ON",
        "-DBUILDNAME:STRING=\"%s\"" % (BUILD_NAME),
        "-DTETGEN_DIR:PATH=%s" % (TETGEN),
        "-DTETGEN_LIBRARY:FILEPATH=%s/tetgen.lib" % (TETGEN),
        "-DTEEM_DIR:PATH=%s/%s64" % (THIRDPARTY, build_type),
        "%s" % (SRC, )])

try:
    os.unlink("%s/scirun.exe" %(build_type))
except OSError:
    logging.info("I tried to delete the SCIRun executable, but couldn't for some"
                 " reason.  This probably just means you're building SCIRun for"
                 " the first time.  If not, then the file is probably in use.")

# The MS toolchain is pretty unstable, as these things go.  The compiler or
# linker may crash arbitrarily while building SCIRun.  Usually a second run of
# the build gets by the error.  So we try a couple builds before we give up, in
# hopes that it's a system error.
# If we see two build fails in a row, it's probably a developer error, so we
# bail out.
try: 
    build("SCIRUN_CORE.sln", "ALL_BUILD")
except subprocess.CalledProcessError:
    print "First build failed .. trying one more before we give up."
    try:
        msbuild("SCIRUN_CORE.sln", "ALL_BUILD")
    except subprocess.CalledProcessError:
        print "Two builds failed in a row .. bailing out."
        sys.exit(1)


# turn this into a check_call to build the installer.
#C:\Program Files (x86)\Inno Setup 5\iscc.exe" %SRC%/InstallUtils/windows/scirun64.iss

# Copy some files over so we can run SCIRun straight from the build directory.
subprocess.call(["robocopy", "%s/%s/bin/" % (THIRDPARTY, BUILDSTR),
                 "%s" % (build_type),
                 "/S", "/NP", "/NJS", "/NJH"], shell=True)
subprocess.call(["robocopy", "%s/%s/lib/" % (THIRDPARTY, BUILDSTR),
                 "lib/", "/S", "/NP", "/NJS", "/NJH"], shell=True)

# Can't seem to get the `copy' command to do what we want ...
for fn in os.listdir("lib/%s/" % (build_type)):
    if fn.endswith(".dll"):
        print "Copying", fn, "..."
        shutil.copy("lib/%s/%s" % (build_type, fn), build_type)

fp = open("ctest-log.txt", "w")
if options.run_tests:
    subprocess.check_call(["ctest", "-V", "-D", "Experimental",
                           "-A", "CMakeCache.txt",
                           "--build-nocmake", "--build-noclean"],
                           stderr=subprocess.STDOUT, stdout=fp
                         )
fp.close()
