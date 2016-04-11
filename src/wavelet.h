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


#ifndef wavelet_h_
#define wavelet_h_

#include "cArrayOp.h"

#include <vector>

#include "itkImage.h"

#include "boost/shared_ptr.hpp"

namespace gth818n
{
  /*============================================================
    wavelet transform */
  template< typename TDouble >
  typename cArray2D< TDouble >::PointerType waveletTransform(typename cArray2D< TDouble >::PointerType array, \
                                                             bool forwardTrans, \
                                                             int daubechiesType = 4, \
                                                             char transformType = 0);

  template< typename TDouble >
  typename cArray3D< TDouble >::PointerType waveletTransform(typename cArray3D< TDouble >::PointerType array, \
                                                             bool forwardTrans, \
                                                             int daubechiesType = 4, \
                                                             char transformType = 0);



  /*==================================================
    DWT2 */
  template< typename TDouble >
  typename cArray2D< TDouble >::PointerType dwt2(typename cArray2D< TDouble >::PointerType array);

  /*==================================================
    IDWT2 */
  template< typename TDouble >
  typename cArray2D< TDouble >::PointerType idwt2(typename cArray2D< TDouble >::PointerType dwtCoef);


 /*==================================================
    DWT3 */
  template< typename TDouble >
  typename cArray3D< TDouble >::PointerType dwt3(typename cArray3D< TDouble >::PointerType array);

  /*==================================================
    IDWT3 */
  template< typename TDouble >
  typename cArray3D< TDouble >::PointerType idwt3(typename cArray3D< TDouble >::PointerType dwtCoef);


  /*==================================================
    DWT2
    For Itk image type
  */
  template< typename TDouble >
  typename itk::Image<TDouble, 2>::Pointer dwt2(typename itk::Image<TDouble, 2>::Pointer array);

  /*==================================================
    IDWT2
    For Itk image type
  */
  template< typename TDouble >
  typename itk::Image<TDouble, 2>::Pointer idwt2(typename itk::Image<TDouble, 2>::Pointer dwtCoef);


 /*==================================================
    DWT3
    For Itk image type
 */
  template< typename TDouble >
  typename itk::Image<TDouble, 3>::Pointer dwt3(typename itk::Image<TDouble, 3>::Pointer array);

  /*==================================================
    IDWT3
    For Itk image type
  */
  template< typename TDouble >
  typename itk::Image<TDouble, 3>::Pointer idwt3(typename itk::Image<TDouble, 3>::Pointer dwtCoef);



  /*==================================================
    make the size right for DWT*/
  template< typename TPixel >
  typename cArray2D< TPixel >::PointerType padToRightSize(typename cArray2D< TPixel >::PointerType array, TPixel padValue = 0);

  template< typename TPixel >
  typename cArray3D< TPixel >::PointerType padToRightSize(typename cArray3D< TPixel >::PointerType array, TPixel padValue = 0);




  /*==================================================
    hard threshold by magnitude:
    all the array values, whose abs is less than t, is set to 0.0 */
  template< typename TDouble >
  void hardThresholdByMagnitude(typename cArray2D< TDouble >::PointerType array, double t);

  template< typename TDouble >
  void hardThresholdByMagnitude(typename cArray3D< TDouble >::PointerType array, double t);



  /*==================================================
    hard threshold by percentage, set r% coef to 0 */
  template< typename TDouble >
  void hardThresholdByPercentage(typename cArray2D< TDouble >::PointerType array, double r = 0.0);

  template< typename TDouble >
  void hardThresholdByPercentage(typename cArray3D< TDouble >::PointerType array, double r = 0.0);


  /**
   * 2D
   */
  template< typename TDouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  coefArrayToBandStructure(const typename cArray2D< TDouble >::PointerType coef,  \
                           char transformType = 0);

  /**
   * 2D, copy of above, just change name
   */
  template< typename TDouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  wArrayToBandStructure(const typename cArray2D< TDouble >::PointerType coef, \
                           char transformType = 0);

  /**
   * 2D, copy of above, for itkImage
   */
  template< typename TDouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  wImageToBandStructure(const typename itk::Image< TDouble, 2 >::Pointer wimage, \
                        char transformType = 0);

  /**
   * 3D
   */
  template< typename TDouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  coefArrayToBandStructure(const typename cArray3D< TDouble >::PointerType coef, \
                           char transformType = 0);

  /**
   * 3D, copy of above, just change name
   */
  template< typename TDouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  wArrayToBandStructure(const typename cArray3D< TDouble >::PointerType coef, \
                           char transformType = 0);

  /**
   * 3D, copy of above, for iktImage
   */
  template< typename TDouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  wImageToBandStructure(const typename itk::Image< TDouble, 3 >::Pointer wimage, \
                        char transformType = 0);


  template< typename TDouble >
  typename cArray2D< TDouble >::PointerType
  bandStructureToCoefArray2D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                              long nx, long ny,                         \
                              char transformType = 0);

  // copy of above, just change name
  template< typename TDouble >
  typename cArray2D< TDouble >::PointerType
  bandStructureToWArray2D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                           long nx, long ny,                            \
                           char transformType = 0);


  template< typename TDouble >
  typename cArray3D< TDouble >::PointerType
  bandStructureToCoefArray3D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                              long nx, long ny, long nz,                \
                              char transformType = 0);

  // copy of above, just change name
  template< typename TDouble >
  typename cArray3D< TDouble >::PointerType
  bandStructureToWArray3D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                           long nx, long ny, long nz,                   \
                           char transformType = 0);

  /**
   * copy of above, for itk image
   */
  template< typename TDouble >
  typename itk::Image< TDouble, 3 >::Pointer
  bandStructureToWimage( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                         long nx, long ny, long nz,                     \
                         char transformType = 0);

} // namespace



#include "wavelet.hxx"

#endif
