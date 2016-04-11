/*******************************************************************************/
/*                                                                             */
/*  Copyright (c) 2010, Yi Gao                                                 */
/*  gaoyi@gatech.edu                                                           */
/*                                                                             */
/*                                                                             */
/*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
/*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   */
/*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    */
/*  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
/*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    */
/*  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        */
/*  DEALINGS IN THE SOFTWARE.                                                  */
/*                                                                             */
/*  See the README.md and COPYING files for usage and copyright information.   */
/*                                                                             */
/*******************************************************************************/

#ifndef cArrayOp_h_
#define cArrayOp_h_

#include "itkImage.h"

#include <vnl/vnl_vector.h>


#include "basicHeaders.h"


#include "cArray2D.h"
#include "cArray3D.h"


namespace gth818n
{
  /* IO with images   */
  template< typename T >
  typename cArray2D< T >::PointerType readImageToArray2(const char *fileName);

  template< typename T >
  typename cArray3D< T >::PointerType readImageToArray3(const char *fileName);

  template< typename T >
  void saveAsImage2(typename cArray2D< T >::PointerType array2, const char *fileName);

  template< typename T >
  void saveAsImage3(typename cArray3D< T >::PointerType array3, const char *fileName);


  /* IO: image file <===> itk::Image  */
  template< typename T > typename itk::Image< T, 2 >::Pointer readImage2(const char *fileName);
  template< typename T > typename itk::Image< T, 3 >::Pointer readImage3(const char *fileName);
  template< typename T, unsigned int dim > typename itk::Image< T, dim >::Pointer readImage(const char *fileName);

  template< typename T > void writeImage2(typename itk::Image< T, 2 >::Pointer img, const char *fileName);
  template< typename T > void writeImage3(typename itk::Image< T, 3 >::Pointer img, const char *fileName);
  template< typename T, unsigned int dim > void writeImage(typename itk::Image< T, dim >::Pointer img, const char *fileName);


  /*  IO: cArray <===> itk::Image  */
  template< typename T >
  typename cArray2D< T >::PointerType itkImageToArray2(typename itk::Image< T, 2 >::Pointer image);

  template< typename T >
  typename cArray3D< T >::PointerType itkImageToArray3(typename itk::Image< T, 3 >::Pointer image);

  template< typename T >
  typename itk::Image< T, 2 >::Pointer cArray2ToItkImage(typename cArray2D< T >::PointerType array2);

  template< typename T >
  typename itk::Image< T, 3 >::Pointer cArray3ToItkImage(typename cArray3D< T >::PointerType array3);


  /* crop array  */
  template< typename T >
  typename cArray2D< T >::PointerType cropArray(typename cArray2D< T >::PointerType array2, \
                                                long xMin, long xMax, long yMin, long yMax);

  template< typename T >
  typename cArray3D< T >::PointerType cropArray(typename cArray3D< T >::PointerType array3, \
                                                long xMin, long xMax, long yMin, long yMax, long zMin, long zMax);


  /* const pad array  */
  template< typename T >
  typename cArray2D< T >::PointerType constPadArray(typename cArray2D< T >::PointerType array2, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum);

  template< typename T >
  typename cArray3D< T >::PointerType constPadArray(typename cArray3D< T >::PointerType array3, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum, \
                                                    long zHeadNum, long zTailNum);
}

#include "cArrayOp.hxx"

#endif
/* definition in .hxx*/
