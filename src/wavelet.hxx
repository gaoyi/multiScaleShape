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


#ifndef wavelet_txx_
#define wavelet_txx_


#include "wavelet.h"

#include <csignal>

#include "gsl/gsl_wavelet.h"

#include <algorithm>



namespace gth818n
{
  static void dwt_step(const gsl_wavelet* w, double* a, size_t stride, size_t n, gsl_wavelet_direction dir, gsl_wavelet_workspace* work)
  {
    // copy from gsl source code
#define ELEMENT(a,stride,i) ((a)[(stride)*(i)])

    double ai, ai1;
    size_t i, ii;
    size_t jf;
    size_t k;
    size_t n1, ni, nh, nmod;

    for (i = 0; i < work->n; i++)
      {
        work->scratch[i] = 0.0;
      }

    nmod = w->nc * n;
    nmod -= w->offset;            /* center support */

    n1 = n - 1;
    nh = n >> 1;

    if (dir == gsl_wavelet_forward)
      {
        for (ii = 0, i = 0; i < n; i += 2, ii++)
          {
            double h = 0, g = 0;

            ni = i + nmod;

            for (k = 0; k < w->nc; k++)
              {
                jf = n1 & (ni + k);
                h += w->h1[k] * ELEMENT(a, stride, jf);
                g += w->g1[k] * ELEMENT(a, stride, jf);
              }

            work->scratch[ii] += h;
            work->scratch[ii + nh] += g;
          }
      }
    else
      {
        for (ii = 0, i = 0; i < n; i += 2, ii++)
          {
            ai = ELEMENT(a, stride, ii);
            ai1 = ELEMENT(a, stride, ii + nh);
            ni = i + nmod;
            for (k = 0; k < w->nc; k++)
              {
                jf = (n1 & (ni + k));
                work->scratch[jf] += (w->h2[k] * ai + w->g2[k] * ai1);
              }
          }
      }

    for (i = 0; i < n; i++)
      {
        ELEMENT(a, stride, i) = work->scratch[i];
      }
  }


  /*============================================================
    wavelet transform */
  template< typename Tdouble >
  typename cArray2D< Tdouble >::PointerType waveletTransform(typename cArray2D< Tdouble >::PointerType array, \
                                                             bool forwardTrans, \
                                                             int daubechiesType, \
                                                             char transformType)
  {
    long nx = array->getSizeX();
    double logn = log(nx)/log(2.0);
    if ( logn - floor(logn) > 1e-6 )
      {
        std::cerr<<"x-dimension is not exp of 2. abort\n";
        raise(SIGABRT);
      }

    long ny = array->getSizeY();

    logn = log(ny)/log(2.0);
    if ( logn - floor(logn) > 1e-6 )
      {
        std::cerr<<"y-dimension is not exp of 2. abort\n";
        raise(SIGABRT);
      }

    long n = nx<ny?nx:ny;


    typename cArray2D< Tdouble >::PointerType dwtCoef(new cArray2D< Tdouble >(*array) );
    if (n < 2)
      {
        return dwtCoef;
      }


    // cArray2D is x-first y-second storing structure

    // choose wavelet
    gsl_wavelet* w = gsl_wavelet_alloc(gsl_wavelet_daubechies, daubechiesType);


    double* dataPointer = dwtCoef->getDataPointer();

    // performing 1D DWT for each *row*, *INPLACE* operation
    gsl_wavelet_workspace* workX = gsl_wavelet_workspace_alloc(nx);
    gsl_wavelet_workspace* workY = gsl_wavelet_workspace_alloc(ny);

    if (0 == transformType)
      {
        /* 0 == transformType:
           1 step in x dir, then 1 step in y dir
           This is called the "NON-standard transform in gsl reference. */
        if ( forwardTrans )
          {
            gsl_wavelet_direction dir = gsl_wavelet_forward;

            long nxx = nx;
            long nyy = ny;

            for (long i = n; i >= 2; i >>= 1)
              {
                for (long iy = 0; iy < nyy; ++iy)       /* for every row iy */
                  {
                    double* rowPtr = iy*nx + dataPointer;
                    dwt_step(w, rowPtr, 1, nxx, dir, workX);
                  }

                for (long ix = 0; ix < nxx; ++ix)       /* for every column ix */
                  {
                    double* colPtr = ix + dataPointer;
                    dwt_step(w, colPtr, nx, nyy, dir, workY);
                  }

                nxx >>= 1;
                nyy >>= 1;
              }
          }
        else
          {
            gsl_wavelet_direction dir = gsl_wavelet_backward;

            /* Compute the size of smallest decomposition band, which
             contains the largest scale coef and the next largest
             wavelet coef. ( the one on the top left corner).

            If the data is square, then it's 2x2, but in case of
            non-square, it's a rectancgle with one side is 2 and the
            other 4, 8 or larger. */

            long nxx = nx;
            long nyy = ny;
            for (long i = n; i > 2; i >>= 1) // > 2, NOT >= 2
              {
                nxx >>= 1;
                nyy >>= 1;
              }

            for (long i = 2; i <= n; i <<= 1)
              {
                for (long ix = 0; ix < nxx; ++ix)       /* for every column ix */
                  {
                    double* colPtr = ix + dataPointer;
                    dwt_step(w, colPtr, nx, nyy, dir, workY);
                  }
                for (long iy = 0; iy < nyy; ++iy)       /* for every row iy */
                  {
                    double* rowPtr = iy*nx + dataPointer;
                    dwt_step(w, rowPtr, 1, nxx, dir, workX);
                  }

                nxx <<= 1;
                nyy <<= 1;
              }
          }
      }
    else if (1 == transformType)
      {
        /* 1 == transformType:
           all steps in x dir, then all steps in y dir
           This is called the "standard transform in gsl reference. */

        if ( forwardTrans )
          {
            for (long iy = 0; iy <= ny-1; ++iy)
              {
                double* rowPtr = iy*nx + dataPointer;
                gsl_wavelet_transform_forward(w, rowPtr, 1, nx, workX);
              }

            for (long ix = 0; ix <= nx-1; ++ix)
              {
                double* colPtr = ix + dataPointer;
                gsl_wavelet_transform_forward(w, colPtr, nx, ny, workY);
              }
          }
        else
          {
            for (long iy = 0; iy <= ny-1; ++iy)
              {
                double* rowPtr = iy*nx + dataPointer;
                gsl_wavelet_transform_inverse(w, rowPtr, 1, nx, workX);
              }

            for (long ix = 0; ix <= nx-1; ++ix)
              {
                double* colPtr = ix + dataPointer;
                gsl_wavelet_transform_inverse(w, colPtr, nx, ny, workY);
              }
          }
      }
    else
      {
        std::cerr<<"transformType should be 0, or 1. abort.\n";
        raise(SIGABRT);
      }



    gsl_wavelet_free(w);
    gsl_wavelet_workspace_free(workX);
    gsl_wavelet_workspace_free(workY);

    return dwtCoef;
  }




  /*============================================================*/
  template< typename Tdouble >
  typename cArray3D< Tdouble >::PointerType waveletTransform(typename cArray3D< Tdouble >::PointerType array, \
                                                             bool forwardTrans, \
                                                             int daubechiesType, \
                                                             char transformType)
  {
    long nx = array->getSizeX();
    double logn = log(nx)/log(2.0);
    if ( logn - floor(logn) > 1e-6 )
      {
        std::cerr<<"x-dimension is not exp of 2. abort\n";
        raise(SIGABRT);
      }

    long ny = array->getSizeY();
    logn = log(ny)/log(2.0);
    if ( logn - floor(logn) > 1e-6 )
      {
        std::cerr<<"y-dimension is not exp of 2. abort\n";
        raise(SIGABRT);
      }

    long nz = array->getSizeZ();
    logn = log(nz)/log(2.0);
    if ( logn - floor(logn) > 1e-6 )
      {
        std::cerr<<"z-dimension is not exp of 2. abort\n";
        raise(SIGABRT);
      }


    long n = nx<ny?nx:ny;
    n = n<nz?n:nz;


    // *INPLACE* operation, thus dwtCoef will store all the coefficients at the end.
    typename cArray3D< Tdouble >::PointerType dwtCoef(new cArray3D< Tdouble >(*array) );
    if (n < 2)
      {
        return dwtCoef;
      }


    // cArray3D is x-first y-second z-second storing structure

    // choose wavelet
    gsl_wavelet* w = gsl_wavelet_alloc(gsl_wavelet_daubechies, daubechiesType);


    double* dataPointer = dwtCoef->getDataPointer();

    // performing 1D DWT for each *row*, *INPLACE* operation
    gsl_wavelet_workspace* workX = gsl_wavelet_workspace_alloc(nx);
    gsl_wavelet_workspace* workY = gsl_wavelet_workspace_alloc(ny);
    gsl_wavelet_workspace* workZ = gsl_wavelet_workspace_alloc(nz);

    long nxy = nx*ny;

    if (0 == transformType)
      {
        /* 0 == transformType:
           1 step in x dir, then 1 step in y dir
           This is called the "NON-standard transform in gsl reference. */

        if ( forwardTrans )
          {
            gsl_wavelet_direction dir = gsl_wavelet_forward;

            long nxx = nx;
            long nyy = ny;
            long nzz = nz;

            for (long i = n; i >= 2; i >>= 1)
              {
                for (long iy = 0; iy < nyy; ++iy)
                  {
                    for (long iz = 0; iz < nzz; ++iz)
                      {
                        double* rowPtr = iz*nxy + iy*nx + dataPointer;
                        dwt_step(w, rowPtr, 1, nxx, dir, workX);
                      }
                  }


                for (long ix = 0; ix < nxx; ++ix)
                  {
                    for (long iz = 0; iz < nzz; ++iz)
                      {
                        double *colPtr = iz*nxy + ix + dataPointer;
                        dwt_step(w, colPtr, nx, nyy, dir, workY);
                      }
                  }



                for (long ix = 0; ix < nxx; ++ix)
                  {
                    for (long iy = 0; iy < nyy; ++iy)
                      {
                        double *zPtr = iy*nx + ix + dataPointer;
                        dwt_step(w, zPtr, nxy, nzz, dir, workZ);
                      }
                  }

                nxx >>= 1;
                nyy >>= 1;
                nzz >>= 1;
              }
          }
        else
          {
            gsl_wavelet_direction dir = gsl_wavelet_backward;

            /* Compute the size of smallest decomposition band, which
             contains the largest scale coef and the next largest
             wavelet coef. ( the one on the top left corner).

            If the data is square, then it's 2x2x2, but in case of
            non-square, it's a rectancgle with one side is 2 and the
            other 4, 8 or larger. */

            long nxx = nx;
            long nyy = ny;
            long nzz = nz;
            for (long i = n; i > 2; i >>= 1) // > 2, NOT >= 2
              {
                nxx >>= 1;
                nyy >>= 1;
                nzz >>= 1;
              }


            for (long i = 2; i <= n; i <<= 1)
              {
                // z first, since in forward trans, z last
                for (long ix = 0; ix < nxx; ++ix)
                  {
                    for (long iy = 0; iy < nyy; ++iy)
                      {
                        double *zPtr = iy*nx + ix + dataPointer;
                        dwt_step(w, zPtr, nxy, nzz, dir, workZ);
                      }
                  }

                // column 2nd, since in forward trans, column 2nd
                for (long ix = 0; ix < nxx; ++ix)
                  {
                    for (long iz = 0; iz < nzz; ++iz)
                      {
                        double *colPtr = iz*nxy + ix + dataPointer;
                        dwt_step(w, colPtr, nx, nyy, dir, workY);
                      }
                  }

                // row 3rd, since in forward trans, row 1st
                for (long iy = 0; iy < nyy; ++iy)
                  {
                    for (long iz = 0; iz < nzz; ++iz)
                      {
                        double* rowPtr = iz*nxy + iy*nx + dataPointer;
                        dwt_step(w, rowPtr, 1, nxx, dir, workX);
                      }
                  }


                nxx <<= 1;
                nyy <<= 1;
                nzz <<= 1;
              }
          }
      }
    else if ( 1 == transformType )
      {
        /* 1 == transformType:
           all steps in x dir, then all steps in y dir
           This is called the "standard transform in gsl reference. */

        if (forwardTrans)
          {
            // forward transform
            for (long iy = 0; iy <= ny-1; ++iy)
              {
                for (long iz = 0; iz <= nz-1; ++iz)
                  {
                    double* rowPtr = iz*nxy + iy*nx + dataPointer;
                    gsl_wavelet_transform_forward(w, rowPtr, 1, nx, workX);
                  }
              }

            for (long ix = 0; ix <= nx-1; ++ix)
              {
                for (long iz = 0; iz <= nz-1; ++iz)
                  {
                    double *colPtr = iz*nxy + ix + dataPointer;
                    gsl_wavelet_transform_forward(w, colPtr, nx, ny, workY);
                  }
              }

            for (long ix = 0; ix <= nx-1; ++ix)
              {
                for (long iy = 0; iy <= ny-1; ++iy)
                  {
                    double *zPtr = iy*nx + ix + dataPointer;
                    gsl_wavelet_transform_forward(w, zPtr, nxy, ny, workZ);
                  }
              }
          }
        else
          {
            // back transform

            for (long iy = 0; iy <= ny-1; ++iy)
              {
                for (long iz = 0; iz <= nz-1; ++iz)
                  {
                    double* rowPtr = iz*nxy + iy*nx + dataPointer;
                    gsl_wavelet_transform_inverse(w, rowPtr, 1, nx, workX);
                  }
              }

            for (long ix = 0; ix <= nx-1; ++ix)
              {
                for (long iz = 0; iz <= nz-1; ++iz)
                  {
                    double *colPtr = iz*nxy + ix + dataPointer;
                    gsl_wavelet_transform_inverse(w, colPtr, nx, ny, workY);
                  }
              }

            for (long ix = 0; ix <= nx-1; ++ix)
              {
                for (long iy = 0; iy <= ny-1; ++iy)
                  {
                    double *zPtr = iy*nx + ix + dataPointer;
                    gsl_wavelet_transform_inverse(w, zPtr, nxy, ny, workZ);
                  }
              }
          }

      }



    gsl_wavelet_free(w);
    gsl_wavelet_workspace_free(workX);
    gsl_wavelet_workspace_free(workY);
    gsl_wavelet_workspace_free(workZ);

    return dwtCoef;
  }



  /*==================================================
    make the size right for DWT, 2D */
  template< typename TPixel >
  typename cArray2D< TPixel >::PointerType padToRightSize(typename cArray2D< TPixel >::PointerType array, TPixel padValue)
  {
    double logn = 0;

    // find right nx
    long nx = array->getSizeX();
    logn = log(nx)/log(2.0);
    long newNx = 0;
    if ( logn - floor(logn) < 1e-6 )
      {
        newNx = nx;
      }
    else
      {
        newNx = (long)pow(2.0, floor(log(nx)/log(2)) + 1.0);
      }

    // find right ny
    long ny = array->getSizeY();
    logn = log(ny)/log(2.0);
    long newNy = 0;
    if ( logn - floor(logn) < 1e-6 )
      {
        newNy = ny;
      }
    else
      {
        newNy = (long)pow(2.0, floor(log(ny)/log(2)) + 1.0);
      }

    long xHeadNum = (long)floor((newNx - nx)/2.0);
    long xTailNum = newNx - nx - xHeadNum;

    long yHeadNum = (long)floor((newNy - ny)/2.0);
    long yTailNum = newNy - ny - yHeadNum;

    typename cArray2D< TPixel >::PointerType paddedArray = constPadArray<TPixel>(array, padValue, xHeadNum, xTailNum, yHeadNum, yTailNum);


    return paddedArray;
  }


  /*==================================================
    make the size right for DWT, 3D */
  template< typename TPixel >
  typename cArray3D< TPixel >::PointerType padToRightSize(typename cArray3D< TPixel >::PointerType array, TPixel padValue)
  {
    double logn = 0;

    // find right nx
    long nx = array->getSizeX();
    logn = log(nx)/log(2.0);
    long newNx = 0;
    if ( logn - floor(logn) < 1e-6 )
      {
        newNx = nx;
      }
    else
      {
        newNx = (long)pow(2.0, floor(log(nx)/log(2)) + 1.0);
      }

    // find right ny
    long ny = array->getSizeY();
    logn = log(ny)/log(2.0);
    long newNy = 0;
    if ( logn - floor(logn) < 1e-6 )
      {
        newNy = ny;
      }
    else
      {
        newNy = (long)pow(2.0, floor(log(ny)/log(2)) + 1.0);
      }

    // find right nz
    long nz = array->getSizeZ();
    logn = log(nz)/log(2.0);
    long newNz = 0;
    if ( logn - floor(logn) < 1e-6 )
      {
        newNz = nz;
      }
    else
      {
        newNz = (long)pow(2.0, floor(log(nz)/log(2)) + 1.0);
      }


    long xHeadNum = (long)floor((newNx - nx)/2.0);
    long xTailNum = newNx - nx - xHeadNum;

    long yHeadNum = (long)floor((newNy - ny)/2.0);
    long yTailNum = newNy - ny - yHeadNum;

    long zHeadNum = (long)floor((newNz - nz)/2.0);
    long zTailNum = newNz - nz - zHeadNum;

    typename cArray3D< TPixel >::PointerType paddedArray = constPadArray<TPixel>(array, padValue, \
                                                                                 xHeadNum, xTailNum, yHeadNum, yTailNum, zHeadNum, zTailNum);


    return paddedArray;
  }




  /*==================================================
    DWT2
    For Itk image type
  */
  template< typename TDouble >
  typename itk::Image<TDouble, 2>::Pointer dwt2(typename itk::Image<TDouble, 2>::Pointer array)
  {
    /*
      1. itk image to carray
      2. pad carray
      3. dwt of padded carray
      4. coef carray to itk image
      5. copy information from array so that it's not lost(though those info does not fit wavelet coefs)
    */

    typedef cArray2D< TDouble > TArray;
    typedef itk::Image< TDouble, 2 > TImage;

    // 1.
    typename TArray::PointerType a = itkImageToArray2<TDouble>(array);

    // 2.
    typename TArray::PointerType aRightSize = padToRightSize<TDouble>(a);

    // 3.
    a.reset();// free memory of im by dropping its ref count by 1
    typename TArray::PointerType dwtCoef = dwt2< TDouble >(aRightSize);

    // 4.
    typename TImage::Pointer dwtCoefImg = cArray2ToItkImage<TDouble>(dwtCoef);

    // 5.
    dwtCoefImg->SetSpacing(array->GetSpacing());
    dwtCoefImg->CopyInformation(array);

    return dwtCoefImg;
  }

  /*==================================================
    IDWT2
    For Itk image type
  */
  template< typename TDouble >
  typename itk::Image<TDouble, 2>::Pointer idwt2(typename itk::Image<TDouble, 2>::Pointer dwtCoef)
  {
    /*
      Assume dwtCoef has the right size!!! (Otherwise where it is from???)

      1. itk image to carray
      2. idwt of carray
      3. carray to itk image
    */

    typedef cArray2D< TDouble > TArray;
    typedef itk::Image< TDouble, 2 > TImage;

    // 1.
    typename TArray::PointerType a = itkImageToArray2<TDouble>(dwtCoef);

    // 2.
    typename TArray::PointerType idwtRslt = idwt2< TDouble >(a);

    // 3.
    typename TImage::Pointer img = cArray2ToItkImage<TDouble>(idwtRslt);

    // 4.
    img->SetSpacing(dwtCoef->GetSpacing());
    img->CopyInformation(dwtCoef);


    return img;
  }


 /*==================================================
    DWT3
    For Itk image type
 */
  template< typename TDouble >
  typename itk::Image<TDouble, 3>::Pointer dwt3(typename itk::Image<TDouble, 3>::Pointer array)
  {
    /*
      1. itk image to carray
      2. pad carray
      3. dwt of padded carray
      4. coef carray to itk image
    */

    typedef cArray3D< TDouble > TArray;
    typedef itk::Image< TDouble, 3 > TImage;

    // 1.
    typename TArray::PointerType a = itkImageToArray3<TDouble>(array);

    // 2.
    typename TArray::PointerType aRightSize = padToRightSize<TDouble>(a);

    // 3.
    a.reset();// free memory of im by dropping its ref count by 1
    typename TArray::PointerType dwtCoef = dwt3< TDouble >(aRightSize);

    // 4.
    typename TImage::Pointer dwtCoefImg = cArray3ToItkImage<TDouble>(dwtCoef);

    //5.
    dwtCoefImg->SetSpacing(array->GetSpacing());
    dwtCoefImg->CopyInformation(array);


    return dwtCoefImg;
  }

  /*==================================================
    IDWT3
    For Itk image type
  */
  template< typename TDouble >
  typename itk::Image<TDouble, 3>::Pointer idwt3(typename itk::Image<TDouble, 3>::Pointer dwtCoef)
  {
    /*
      Assume dwtCoef has the right size!!! (Otherwise where it is from???)

      1. itk image to carray
      2. idwt of carray
      3. carray to itk image
    */

    typedef cArray3D< TDouble > TArray;
    typedef itk::Image< TDouble, 3 > TImage;

    // 1.
    typename TArray::PointerType a = itkImageToArray3<TDouble>(dwtCoef);

    // 2.
    typename TArray::PointerType idwtRslt = idwt3< TDouble >(a);

    // 3.
    typename TImage::Pointer img = cArray3ToItkImage<TDouble>(idwtRslt);

    // 4.
    img->SetSpacing(dwtCoef->GetSpacing());
    img->CopyInformation(dwtCoef);


    return img;
  }



  /*==================================================
    DWT2 */
  template< typename Tdouble >
  typename cArray2D< Tdouble >::PointerType dwt2(typename cArray2D< Tdouble >::PointerType array)
  {
    bool forwardTrans = true;

    return waveletTransform< Tdouble >(array, forwardTrans);
  }


  /*==================================================
    IDWT2 */
  template< typename Tdouble >
  typename cArray2D< Tdouble >::PointerType idwt2(typename cArray2D< Tdouble >::PointerType dwtCoef)
  {
    bool forwardTrans = false;

    return waveletTransform< Tdouble >(dwtCoef, forwardTrans);
  }



  /*==================================================
    DWT3 */
  template< typename Tdouble >
  typename cArray3D< Tdouble >::PointerType dwt3(typename cArray3D< Tdouble >::PointerType array)
  {
    bool forwardTrans = true;

    return waveletTransform< Tdouble >(array, forwardTrans);
  }


  /*==================================================
    IDWT3 */
  template< typename Tdouble >
  typename cArray3D< Tdouble >::PointerType idwt3(typename cArray3D< Tdouble >::PointerType dwtCoef)
  {
    bool forwardTrans = false;

    return waveletTransform< Tdouble >(dwtCoef, forwardTrans);
  }



  /*==================================================
    hard threshold by magnitude:
    all the array values, whose abs is less than t, is set to 0.0 */
  template< typename Tdouble >
  void hardThresholdByMagnitude(typename cArray2D< Tdouble >::PointerType array, double t)
  {
    long nx = array->getSizeX();
    long ny = array->getSizeY();
    long n = nx*ny;

    for (long i = 0; i < n; ++i)
      {
        double v = array->get(i);
        double absv = fabs(v);

        array->set(i, absv>=t?v:0.0 );
      }

  }


  template< typename Tdouble >
  void hardThresholdByMagnitude(typename cArray3D< Tdouble >::PointerType array, double t)
  {
    long nx = array->getSizeX();
    long ny = array->getSizeY();
    long nz = array->getSizeZ();
    long n = nx*ny*nz;

    for (long i = 0; i < n; ++i)
      {
        double v = array->get(i);
        double absv = fabs(v);

        array->set(i, absv>=t?v:0.0);
      }

  }


  /*==================================================
    hard threshold by percentage, set r% coef to 0 */
  template< typename Tdouble >
  void hardThresholdByPercentage(typename cArray2D< Tdouble >::PointerType array, double r)
  {
    if (r < 0 || r > 1)
      {
        std::cerr<<"r should be in [0, 1]. abort\n";
        raise(SIGABRT);
      }

    typedef cArray2D< Tdouble > ArrayType;

    long nx = array->getSizeX();
    long ny = array->getSizeY();
    long n = nx*ny;

    // 1. sort the element by there magnitude
    ArrayType m(*array);
    for (long i = 0; i < n; ++i)
      {
        m[i] = fabs(m[i]);
      }

    double* pm = m.getDataPointer();

    std::sort(pm, pm + n);

    // 2. find the magnitude
    long k = static_cast<long>(floor(r*n));
    double t = m[k];

    // 3. clamp by t
    hardThresholdByMagnitude< Tdouble >(array, t);
  }


  template< typename Tdouble >
  void hardThresholdByPercentage(typename cArray3D< Tdouble >::PointerType array, double r)
  {
    if (r < 0 || r > 1)
      {
        std::cerr<<"r should be in [0, 1]. abort\n";
        raise(SIGABRT);
      }

    typedef cArray3D< Tdouble > ArrayType;

    long nx = array->getSizeX();
    long ny = array->getSizeY();
    long nz = array->getSizeZ();
    long n = nx*ny*nz;

    // 1. sort the element by there magnitude
    ArrayType m(*array);
    for (long i = 0; i < n; ++i)
      {
        m[i] = fabs(m[i]);
      }

    double* pm = m.getDataPointer();

    std::sort(pm, pm + n);

    // 2. find the magnitude
    long k = (long)floor(r*n);
    double t = m[k];

    // 3. clamp by t
    hardThresholdByMagnitude< Tdouble >(array, t);
  }


  /**
   * Arrange the coef into bands.
   * From a 2D/3D matrix, to a std::vector of vnlVectors.
   * Each vnlVector, is the concatinated coef in one band.
   */
  template< typename Tdouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  wArrayToBandStructure(const typename cArray2D< Tdouble >::PointerType coef, \
                        char transformType = 0)
  {
    long nx = coef->getSizeX();
    long ny = coef->getSizeY();

    double log2nX = log(nx)/log(2.0);
    if ( log2nX - floor(log2nX) > 1e-6 )
      {
        std::cerr<<"nx is not power of 2, abort.\n";
        raise(SIGABRT);
      }

    double log2nY = log(ny)/log(2.0);
    if ( log2nY - floor(log2nY) > 1e-6 )
      {
        std::cerr<<"ny is not power of 2, abort.\n";
        raise(SIGABRT);
      }

    typedef vnl_vector< double > VnlDoubleVectorType;

    typedef std::vector< VnlDoubleVectorType > WaveletBandsType;
    typedef boost::shared_ptr< WaveletBandsType > WaveletBandsPointerType;

    WaveletBandsPointerType bands(new WaveletBandsType); // to be returned

    if ( 1 == transformType )
      {
        for (long blockEndX = 1; blockEndX <= nx; blockEndX<<=1)
          {
            for (long blockEndY = 1; blockEndY <= ny; blockEndY<<=1)
              {
                long numOfCoefInThisBand = (blockEndX - (blockEndX>>1))*(blockEndY - (blockEndY>>1));

                VnlDoubleVectorType thisBand(numOfCoefInThisBand );

                long idxInBand = 0;
                for (long inThisBandX = (blockEndX>>1); inThisBandX < blockEndX; ++inThisBandX)
                  {
                    for (long inThisBandY = (blockEndY>>1); inThisBandY < blockEndY; ++inThisBandY)
                      {
                        thisBand[idxInBand] = coef->get(inThisBandX, inThisBandY);
                        ++idxInBand;
                      }
                  }

                bands->push_back(thisBand);
              }
          }
      }
    else if (0 == transformType)
      {
        /* Compute the size of smallest decomposition band, which
           contains the largest scale coef and the next largest
           wavelet coef. ( the one on the top left corner).

           If the data is square, then it's 2x2x2, but in case of
           non-square, it's a rectancgle with one side is 2 and the
           other 4, 8 or larger.

           In the i loop below, it's >= 2, NOT > 2 as in the inverse
           trans. Because we want the coef of the largest scale
           to be a single band, where there the largest and the next
           largest are conbined. */
        long nxx = nx;
        long nyy = ny;

        long n = nx<ny?nx:ny;

        for (long i = n; i >= 2; i >>= 1)
          {
            nxx >>= 1;
            nyy >>= 1;
          }

        //        printf("(nxx, nyy) = (%ld, %ld)\n", nxx, nyy);

        /* the top-left rectangle, containing the coef of the largest
           scale, has the size (nxx, nyy). In case of square, it's (1,
           1). The subsequent process is doubling this square.  */

        // first band
        VnlDoubleVectorType firstBand(nxx*nyy );
        long idxInFirstBand = 0;
        for (long ix = 0; ix < nxx; ++ix)
          {
            for (long iy = 0; iy < nyy; ++iy)
              {
                firstBand[idxInFirstBand] = coef->get(ix, iy);
                ++idxInFirstBand;
              }
          }
        bands->push_back(firstBand);


        nxx <<= 1;
        nyy <<= 1;

        // subsequent bands
        for (long i = 2; i <= n; i <<= 1)
          {
            for (char ia = 1; ia >= 0; --ia )
              {
                for (char ib = 1; ib >= 0; --ib )
                  {
                    if (1 == ia && 1 == ib)
                      {
                        // Correspond to the block of the previous scale (the block to the top AND left)
                        continue;
                      }

                    long blockSizeX = nxx >> 1;
                    long blockSizeY = nyy >> 1;

                    long blockEndX = nxx >> ia;
                    long blockEndY = nyy >> ib;

                    long blockBeginX = blockEndX - blockSizeX;
                    long blockBeginY = blockEndY - blockSizeY;

                    unsigned long numOfCoefInThisBand = blockSizeX*blockSizeY;
                    VnlDoubleVectorType thisBand(numOfCoefInThisBand);

                    long idxInBand = 0;
                    for (long ix = blockBeginX; ix < blockEndX; ++ix)
                      {
                        for (long iy = blockBeginY; iy < blockEndY; ++iy)
                          {
                            thisBand[idxInBand] = coef->get(ix, iy);
                            ++idxInBand;
                          }
                      }

                    bands->push_back(thisBand);
                  }  // ib
              } // ia


            nxx <<= 1;
            nyy <<= 1;
          }
      }
    else
      {
        std::cerr<<"transformType = 1 or 0. abort\n";
        raise(SIGABRT);
      }

    return bands;
  }




  template< typename Tdouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  coefArrayToBandStructure(const typename cArray2D< Tdouble >::PointerType coef, \
                           char transformType)
  {
    return wArrayToBandStructure<Tdouble>(coef, transformType);
  }

  /*============================================================*/
  template< typename Tdouble >
  typename cArray2D< Tdouble >::PointerType
  bandStructureToWArray2D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                           long nx, long ny,                            \
                           char transformType = 0)
  {

    typedef cArray2D< Tdouble > ArrayType;
    typedef vnl_vector< double > VnlDoubleVectorType;

    typename cArray2D< Tdouble >::PointerType array( new ArrayType(nx, ny) );

    if ( 1 == transformType )
      {
        long iBand = 0;
        for (long blockEndX = 1; blockEndX <= nx; blockEndX<<=1)
          {
            for (long blockEndY = 1; blockEndY <= ny; blockEndY<<=1)
              {
                VnlDoubleVectorType thisBand = bands->at(iBand);

                long idxInBand = 0;
                for (long inThisBandX = (blockEndX>>1); inThisBandX < blockEndX; ++inThisBandX)
                  {
                    for (long inThisBandY = (blockEndY>>1); inThisBandY < blockEndY; ++inThisBandY)
                      {
                        array->set(inThisBandX, inThisBandY, thisBand[idxInBand] );
                        ++idxInBand;
                      }
                  }

                ++iBand;
              }
          }
      }
    else if ( 0 == transformType )
      {
        /* Compute the size of smallest decomposition band, which
           contains the largest scale coef and the next largest
           wavelet coef. ( the one on the top left corner).

           If the data is square, then it's 2x2x2, but in case of
           non-square, it's a rectancgle with one side is 2 and the
           other 4, 8 or larger.

           In the i loop below, it's >= 2, NOT > 2 as in the inverse
           trans. Because we want the coef of the largest scale
           to be a single band, where there the largest and the next
           largest are conbined. */
        long nxx = nx;
        long nyy = ny;

        long n = nx<ny?nx:ny;

        unsigned long nb = 0;
        for (long i = n; i >= 2; i >>= 1)
          {
            nxx >>= 1;
            nyy >>= 1;

            ++nb;
          }

        nb = nb*3 + 1;

        if (bands->size() != nb)
          {
            std::cerr<<"bands size wrong. abort\n";
            raise(SIGABRT);
          }

        //        printf("(nxx, nyy) = (%ld, %ld)\n", nxx, nyy);

        /* the top-left rectangle, containing the coef of the largest
           scale, has the size (nxx, nyy). In case of square, it's (1,
           1). The subsequent process is doubling this square.  */

        // first band
        VnlDoubleVectorType firstBand = bands->front();
        if (firstBand.size() != static_cast<unsigned long>(nxx*nyy))
          {
            std::cerr<<"first bands size wrong. abort\n";
            raise(SIGABRT);
          }

        long idxInFirstBand = 0;
        for (long ix = 0; ix < nxx; ++ix)
          {
            for (long iy = 0; iy < nyy; ++iy)
              {
                array->set(ix, iy, firstBand[idxInFirstBand] );
                ++idxInFirstBand;
              }
          }

        nxx <<= 1;
        nyy <<= 1;

        // subsequent bands
        long iband = 0;
        for (long i = 2; i <= n; i <<= 1)
          {
            for (char ia = 1; ia >= 0; --ia )
              {
                for (char ib = 1; ib >= 0; --ib )
                  {
                    if (1 == ia && 1 == ib)
                      {
                        // Correspond to the block of the previous scale (the block to the top AND left)
                        continue;
                      }

                    long blockSizeX = nxx >> 1;
                    long blockSizeY = nyy >> 1;

                    long blockEndX = nxx >> ia;
                    long blockEndY = nyy >> ib;

                    long blockBeginX = blockEndX - blockSizeX;
                    long blockBeginY = blockEndY - blockSizeY;

                    ++iband;
                    VnlDoubleVectorType thisBand = bands->at(iband);
                    if (thisBand.size() != static_cast<unsigned long>(blockSizeX*blockSizeY))
                      {
                        printf("in band %ld, band size is %ld, block size is %uld\n", \
                               i, \
                               static_cast<long>(thisBand.size()), static_cast<unsigned int>(blockSizeX*blockSizeY));

                        std::cerr<<"This bands size wrong. abort\n";
                        raise(SIGABRT);
                      }



                    long idxInBand = 0;
                    for (long ix = blockBeginX; ix < blockEndX; ++ix)
                      {
                        for (long iy = blockBeginY; iy < blockEndY; ++iy)
                          {
                            array->set(ix, iy, thisBand[idxInBand] );
                            ++idxInBand;
                          }
                      }
                  }  // ib
              } // ia


            nxx <<= 1;
            nyy <<= 1;
          }
      }
    else
      {
        std::cerr<<"transformType = 1 or 0. abort\n";
        raise(SIGABRT);
      }



    return array;
  }

  template< typename Tdouble >
  typename cArray2D< Tdouble >::PointerType
  bandStructureToCoefArray2D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                              long nx, long ny, \
                              char transformType)
  {
    return bandStructureToWArray2D<Tdouble>(bands, nx, ny, transformType);
  }


  /*============================================================*/
  template< typename Tdouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  wArrayToBandStructure(const typename cArray3D< Tdouble >::PointerType coef, \
                        char transformType = 0)
  {
    long nx = coef->getSizeX();
    long ny = coef->getSizeY();
    long nz = coef->getSizeZ();

    double log2nX = log(nx)/log(2.0);
    if ( log2nX - floor(log2nX) > 1e-6 )
      {
        std::cerr<<"nx is not power of 2, abort.\n";
        raise(SIGABRT);
      }

    double log2nY = log(ny)/log(2.0);
    if ( log2nY - floor(log2nY) > 1e-6 )
      {
        std::cerr<<"ny is not power of 2, abort.\n";
        raise(SIGABRT);
      }

    double log2nZ = log(nz)/log(2.0);
    if ( log2nZ - floor(log2nZ) > 1e-6 )
      {
        std::cerr<<"nz is not power of 2, abort.\n";
        raise(SIGABRT);
      }


    typedef vnl_vector< double > VnlDoubleVectorType;

    typedef std::vector< VnlDoubleVectorType > WaveletBandsType;
    typedef boost::shared_ptr< WaveletBandsType > WaveletBandsPointerType;


    WaveletBandsPointerType bands(new WaveletBandsType);// to be rtn

    if ( 1 == transformType )
      {
        for (long blockEndX = 1; blockEndX <= nx; blockEndX<<=1)
          {
            for (long blockEndY = 1; blockEndY <= ny; blockEndY<<=1)
              {
                for (long blockEndZ = 1; blockEndZ <= nz; blockEndZ<<=1)
                  {
                    long numOfCoefInThisBand = (blockEndX - (blockEndX>>1))*(blockEndY - (blockEndY>>1))*(blockEndZ - (blockEndZ>>1));

                    VnlDoubleVectorType thisBand( numOfCoefInThisBand );

                    long idxInBand = 0;
                    for (long inThisBandX = (blockEndX>>1); inThisBandX < blockEndX; ++inThisBandX)
                      {
                        for (long inThisBandY = (blockEndY>>1); inThisBandY < blockEndY; ++inThisBandY)
                          {
                            for (long inThisBandZ = (blockEndZ>>1); inThisBandZ < blockEndZ; ++inThisBandZ)
                              {
                                thisBand[idxInBand] = coef->get(inThisBandX, inThisBandY, inThisBandZ);
                                ++idxInBand;
                              }
                          }
                      }

                    bands->push_back(thisBand);
                  }
              }
          }
      }
    else if (0 == transformType)
      {
        /* Compute the size of smallest decomposition band, which
           contains the largest scale coef and the next largest
           wavelet coef. ( the one on the top left corner).

           If the data is square, then it's 2x2x2, but in case of
           non-square, it's a rectancgle with one side is 2 and the
           other 4, 8 or larger.

           In the i loop below, it's >= 2, NOT > 2 as in the inverse
           trans. Because we want the coef of the largest scale
           to be a single band, where there the largest and the next
           largest are conbined. */
        long nxx = nx;
        long nyy = ny;
        long nzz = nz;

        long n = nx<ny?nx:ny;
        n = n<nz?n:nz;

        for (long i = n; i >= 2; i >>= 1)
          {
            nxx >>= 1;
            nyy >>= 1;
            nzz >>= 1;
          }

        //        printf("(nxx, nyy, nzz) = (%ld, %ld, %ld)\n", nxx, nyy, nzz);

        /* the top-left rectangle, containing the coef of the *LARGEST*
         *scale, has the size (nxx, nyy, nzz). In case of square, it's
         *(1, 1, 1). The subsequent process is doubling this
         *square.  */

        // first band
        VnlDoubleVectorType firstBand(nxx*nyy*nzz );
        long idxInFirstBand = 0;
        for (long ix = 0; ix < nxx; ++ix)
          {
            for (long iy = 0; iy < nyy; ++iy)
              {
                for (long iz = 0; iz < nzz; ++iz)
                  {
                    firstBand[idxInFirstBand] = coef->get(ix, iy, iz);
                    ++idxInFirstBand;
                  }
              }
          }
        bands->push_back(firstBand);


        nxx <<= 1;
        nyy <<= 1;
        nzz <<= 1;

        // subsequent bands
        for (long i = 2; i <= n; i <<= 1)
          {
            for (char ia = 1; ia >= 0; --ia )
              {
                for (char ib = 1; ib >= 0; --ib )
                  {
                    for (char ic = 1; ic >= 0; --ic )
                      {
                        if (1 == ia && 1 == ib && 1 == ic)
                          {
                            // Correspond to the block of the previous scale (the block to the top AND left)
                            continue;
                          }

                        long blockSizeX = nxx >> 1;
                        long blockSizeY = nyy >> 1;
                        long blockSizeZ = nzz >> 1;

                        long blockEndX = nxx >> ia;
                        long blockEndY = nyy >> ib;
                        long blockEndZ = nzz >> ic;

                        long blockBeginX = blockEndX - blockSizeX;
                        long blockBeginY = blockEndY - blockSizeY;
                        long blockBeginZ = blockEndZ - blockSizeZ;

                        unsigned long numOfCoefInThisBand = blockSizeX*blockSizeY*blockSizeZ;
                        VnlDoubleVectorType thisBand(numOfCoefInThisBand );

                        long idxInBand = 0;
                        for (long ix = blockBeginX; ix < blockEndX; ++ix)
                          {
                            for (long iy = blockBeginY; iy < blockEndY; ++iy)
                              {
                                for (long iz = blockBeginZ; iz < blockEndZ; ++iz)
                                  {
                                    thisBand[idxInBand] = coef->get(ix, iy, iz);
                                    ++idxInBand;
                                  }
                              }
                          }

                        bands->push_back(thisBand);
                      }  // ic
                  }  // ib
              } // ia


            nxx <<= 1;
            nyy <<= 1;
            nzz <<= 1;
          }
      }
    else
      {
        std::cerr<<"transformType = 1 or 0. abort\n";
        raise(SIGABRT);
      }

    return bands;
  }


  /**
   * 3D, copy of above, for iktImage
   */
  template< typename TDouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  wImageToBandStructure(const typename itk::Image< TDouble, 3 >::Pointer wimage, \
                        char transformType = 0)
  {
    /**
     * 1. itk wImage to wArray
     * 2. wArray to band
     */

    typedef cArray3D< TDouble > TArray;
    typedef itk::Image< TDouble, 3 > itkImage_t;

    // 1.
    typename TArray::Pointer a = itkImageToArray3<TDouble>(wimage);

    // 2.
    return wArrayToBandStructure<TDouble>(a, transformType);
  }



  template< typename Tdouble >
  boost::shared_ptr< std::vector< vnl_vector< double > > >
  coefArrayToBandStructure(const typename cArray3D< Tdouble >::PointerType coef, \
                           char transformType)
  {
    return wArrayToBandStructure<Tdouble>(coef, transformType);
  }


  /*============================================================*/
  template< typename Tdouble >
  typename cArray3D< Tdouble >::PointerType
  bandStructureToWArray3D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                           long nx, long ny, long nz,                   \
                           char transformType = 0)
  {
    typedef cArray3D< Tdouble > ArrayType;
    typedef vnl_vector< double > VnlDoubleVectorType;
    //    typedef boost::shared_ptr< VnlDoubleVectorType > VnlDoubleVectorPointerType;

    typename cArray3D< Tdouble >::PointerType array( new ArrayType(nx, ny, nz) );

    if ( 1 == transformType )
      {
        long iBand = 0;
        for (long blockEndX = 1; blockEndX <= nx; blockEndX<<=1)
          {
            for (long blockEndY = 1; blockEndY <= ny; blockEndY<<=1)
              {
                for (long blockEndZ = 1; blockEndZ <= nz; blockEndZ<<=1)
                  {
                    VnlDoubleVectorType thisBand = bands->at(iBand);

                    long idxInBand = 0;
                    for (long inThisBandX = (blockEndX>>1); inThisBandX < blockEndX; ++inThisBandX)
                      {
                        for (long inThisBandY = (blockEndY>>1); inThisBandY < blockEndY; ++inThisBandY)
                          {
                            for (long inThisBandZ = (blockEndZ>>1); inThisBandZ < blockEndZ; ++inThisBandZ)
                              {
                                array->set(inThisBandX, inThisBandY, inThisBandZ, thisBand[idxInBand] );
                                ++idxInBand;
                              }
                          }
                      }

                    ++iBand;
                  }
              }
          }
      }
    else if ( 0 == transformType )
      {
        /* Compute the size of smallest decomposition band, which
           contains the largest scale coef and the next largest
           wavelet coef. ( the one on the top left corner).

           If the data is square, then it's 2x2x2, but in case of
           non-square, it's a rectancgle with one side is 2 and the
           other 4, 8 or larger.

           In the i loop below, it's >= 2, NOT > 2 as in the inverse
           trans. Because we want the coef of the largest scale
           to be a single band, where there the largest and the next
           largest are conbined. */
        long nxx = nx;
        long nyy = ny;
        long nzz = nz;

        long n = nx<ny?nx:ny;
        n = n<nz?n:nz;

        unsigned long nb = 0;
        for (long i = n; i >= 2; i >>= 1)
          {
            nxx >>= 1;
            nyy >>= 1;
            nzz >>= 1;

            ++nb;
          }

        nb = nb*7 + 1;

        if (bands->size() != nb)
          {
            std::cerr<<"bands size wrong. abort\n";
            raise(SIGABRT);
          }

        /* the top-left rectangle, containing the coef of the largest
           scale, has the size (nxx, nyy, nzz). In case of square, it's (1,
           1, 1). The subsequent process is doubling this square.  */

        // first band
        VnlDoubleVectorType firstBand = bands->front();
        if (firstBand.size() != static_cast<unsigned long>(nxx*nyy*nzz))
          {
            std::cerr<<"first bands size wrong. abort\n";
            raise(SIGABRT);
          }

        long idxInFirstBand = 0;
        for (long ix = 0; ix < nxx; ++ix)
          {
            for (long iy = 0; iy < nyy; ++iy)
              {
                for (long iz = 0; iz < nzz; ++iz)
                  {
                    array->set(ix, iy, iz, firstBand[idxInFirstBand] );
                    ++idxInFirstBand;
                  }
              }
          }

        nxx <<= 1;
        nyy <<= 1;
        nzz <<= 1;

        // subsequent bands
        long iband = 0;
        for (long i = 2; i <= n; i <<= 1)
          {
            for (char ia = 1; ia >= 0; --ia )
              {
                for (char ib = 1; ib >= 0; --ib )
                  {
                    for (char ic = 1; ic >= 0; --ic )
                      {

                        if (1 == ia && 1 == ib && 1 == ic)
                          {
                            // Correspond to the block of the previous scale (the block to the top AND left)
                            continue;
                          }

                        long blockSizeX = nxx >> 1;
                        long blockSizeY = nyy >> 1;
                        long blockSizeZ = nzz >> 1;

                        long blockEndX = nxx >> ia;
                        long blockEndY = nyy >> ib;
                        long blockEndZ = nzz >> ic;

                        long blockBeginX = blockEndX - blockSizeX;
                        long blockBeginY = blockEndY - blockSizeY;
                        long blockBeginZ = blockEndZ - blockSizeZ;

                        ++iband;
                        VnlDoubleVectorType thisBand = bands->at(iband);
                        if (thisBand.size() != static_cast<unsigned long>(blockSizeX*blockSizeY*blockSizeZ))
                          {
                            printf("in band %ld, band size is %ud, block size is %ld\n", \
                                   i, thisBand.size(), blockSizeX*blockSizeY*blockSizeZ);

                            std::cerr<<"This bands size wrong. abort\n";
                            raise(SIGABRT);
                          }


                        long idxInBand = 0;
                        for (long ix = blockBeginX; ix < blockEndX; ++ix)
                          {
                            for (long iy = blockBeginY; iy < blockEndY; ++iy)
                              {
                                for (long iz = blockBeginZ; iz < blockEndZ; ++iz)
                                  {
                                    array->set(ix, iy, iz, thisBand[idxInBand] );
                                    ++idxInBand;
                                  }
                              }
                          }
                      }  // ic
                  }  // ib
              } // ia


            nxx <<= 1;
            nyy <<= 1;
            nzz <<= 1;
          }
      }
    else
      {
        std::cerr<<"transformType = 1 or 0. abort\n";
        raise(SIGABRT);
      }



    return array;
  }

  template< typename TDouble >
  typename cArray3D< TDouble >::PointerType
  bandStructureToCoefArray3D( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                              long nx, long ny, long nz, \
                              char transformType)
  {
    return bandStructureToWArray3D<TDouble>(bands, nx, ny, nz, transformType);
  }


  /**
   * copy of above, for itk image
   */
  template< typename TDouble >
  typename itk::Image< TDouble, 3 >::Pointer
  bandStructureToWimage( const boost::shared_ptr< std::vector< vnl_vector< double > > > bands, \
                         long nx, long ny, long nz,                     \
                         char transformType = 0)
  {
    /**
     * 1. band to wArray
     * 2. wArray to wImage
     */

    // 1.
    typedef cArray3D< TDouble > TArray;
    typename TArray::Pointer a = bandStructureToWArray3D<TDouble>(bands, nx, ny, nz, transformType);

    // 2.
    return cArray3ToItkImage<TDouble>(a);

  }



} // namespace

#endif
