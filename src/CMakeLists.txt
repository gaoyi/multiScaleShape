###############################################################################
#                                                                             #
#  Copyright (c) 2010, Yi Gao                                                 #
#  gaoyi@gatech.edu                                                           #
#                                                                             #
#                                                                             #
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER #
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
#  DEALINGS IN THE SOFTWARE.                                                  #
#                                                                             #
#  See the README.md and COPYING files for usage and copyright information.   #
#                                                                             #
###############################################################################


cmake_minimum_required(VERSION 2.8)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET(projName "Wavelet")
PROJECT(${projName})

find_package (ITK REQUIRED)
include(${ITK_USE_FILE})


# find gsl path
FIND_PATH(GSL_HEADER_PATH include/gsl/gsl_math.h $ENV{GSL_HEADER_PATH} REQUIRED)
IF(NOT GSL_HEADER_PATH)
   MESSAGE(FATAL_ERROR "Please enter the header path of include/gsl/gsl_math.h. the PATH where PATH is set by --prefix when configure gsl")
ENDIF(NOT GSL_HEADER_PATH)

INCLUDE_DIRECTORIES(${GSL_HEADER_PATH}/include/)

# boost
find_path(BOOST_PATH /boost/shared_ptr.hpp $ENV{BOOST_PATH})
if(NOT BOOST_PATH)
   message(FATAL_ERROR "Please enter the path of boost")
endif(NOT BOOST_PATH)

include_directories(${BOOST_PATH})

SET(gslLibs
${GSL_HEADER_PATH}/lib/libgsl.a
)

SET(Libraries
${ITK_LIBRARIES}
${gslLibs}
)



SET(theExe "waveletTest2")
ADD_EXECUTABLE(${theExe} ${theExe}.cxx)
TARGET_LINK_LIBRARIES(${theExe} ${Libraries})

SET(theExe "waveletTest3")
ADD_EXECUTABLE(${theExe} ${theExe}.cxx)
TARGET_LINK_LIBRARIES(${theExe} ${Libraries})

SET(theExe "waveletTestItk3")
ADD_EXECUTABLE(${theExe} ${theExe}.cxx)
TARGET_LINK_LIBRARIES(${theExe} ${Libraries})

SET(theExe "waveletTestItk4")
ADD_EXECUTABLE(${theExe} ${theExe}.cxx)
TARGET_LINK_LIBRARIES(${theExe} ${Libraries})
