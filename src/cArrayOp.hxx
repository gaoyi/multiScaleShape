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


#ifndef cArrayOp_txx_
#define cArrayOp_txx_

#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "cArrayOp.h"

namespace gth818n
{
  /* --------------------------------------------------
     readImage2
  */
  template< typename T > typename itk::Image< T, 2 >::Pointer readImage2(const char *fileName)
  {
    return readImage<T, 2>(fileName);
  }

  /* --------------------------------------------------
     readImage3
  */
  template< typename T > typename itk::Image< T, 3 >::Pointer readImage3(const char *fileName)
  {
    return readImage<T, 3>(fileName);
  }


  /* --------------------------------------------------
     readImage
  */
  template< typename T, unsigned int dim > typename itk::Image< T, dim >::Pointer readImage(const char *fileName)
  {
    typedef itk::Image< T, dim > ImageType;
    typedef itk::ImageFileReader< ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename ImageType::Pointer image;

    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl;
        std::cerr<< err << std::endl;
        raise(SIGABRT);
      }

    return image;
  }


  /*--------------------------------------------------
    writeImage2
  */
  template< typename T > void writeImage2(typename itk::Image< T, 2 >::Pointer img, const char *fileName)
  {
    writeImage<T, 2>(img, fileName);
  }


  /*--------------------------------------------------
    writeImage3
  */
  template< typename T > void writeImage3(typename itk::Image< T, 3 >::Pointer img, const char *fileName)
  {
    writeImage<T, 3>(img, fileName);
  }


  /*--------------------------------------------------
    writeImage
  */
  template< typename T, unsigned int dim > void writeImage(typename itk::Image< T, dim >::Pointer img, const char *fileName)
  {
    typedef itk::Image< T, dim > ImageType;
    typedef itk::ImageFileWriter< ImageType > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(img);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        raise(SIGABRT);
      }
  }


  /*============================================================
    readImageToArray2  */
  template< typename T >
  typename cArray2D< T >::PointerType readImageToArray2(const char *fileName)
  {
    typedef itk::Image< T, 2 > ImageType;
    typedef itk::ImageFileReader< ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename ImageType::Pointer image;

    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl;
        std::cerr<< err << std::endl;
        raise(SIGABRT);
      }

    return itkImageToArray2< T >(image);
  }



  /*============================================================
    readImageToArray3  */
  template< typename T >
  typename cArray3D< T >::PointerType readImageToArray3(const char *fileName)
  {
    typedef itk::Image< T, 3 > ImageType;
    typedef itk::ImageFileReader< ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename ImageType::Pointer image;

    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl;
        std::cerr<< err << std::endl;
        raise(SIGABRT);
      }


    return itkImageToArray3< T >(image);
  }





  /*============================================================
    saveAsImage2  */
  template< typename T >
  void saveAsImage2(typename cArray2D< T >::PointerType array2, const char *fileName)
  {
    typedef itk::Image< T, 2 > ImageType; // 2 means 2D

    typename ImageType::Pointer image = cArray2ToItkImage< T >(array2);

    typedef itk::ImageFileWriter< ImageType > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(image);

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        raise(SIGABRT);
      }

    return;
  }


  /*============================================================
    saveAsImage3  */
  template< typename T >
  void saveAsImage3(typename cArray3D< T >::PointerType array3, const char *fileName)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D

    typename ImageType::Pointer image = cArray3ToItkImage< T >(array3);

    typedef itk::ImageFileWriter< ImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(image);

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        raise(SIGABRT);
      }

    return;
  }

  /*============================================================
    saveAsImage3  */
  template< typename T >
  void saveAsImage3(const cArray3D< T >* array3, const char *fileName)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D

    typename ImageType::Pointer image = cArray3ToItkImage< T >(array3);

    typedef itk::ImageFileWriter< ImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(image);

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        raise(SIGABRT);
      }

    return;
  }



  /* ============================================================
     itkImageToArray2  */
  template< typename T >
  typename cArray2D< T >::PointerType itkImageToArray2(typename itk::Image< T, 2 >::Pointer image)
  {
    typedef itk::Image< T, 2 > ImageType; // 2 means 2D

    typename ImageType::RegionType::SizeType dataSz = (image->GetLargestPossibleRegion()).GetSize();

    typename cArray2D<T>::PointerType array2( new cArray2D< T >(dataSz[0], dataSz[1]) );

    typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;
    ImageIteratorType imIter(image, image->GetLargestPossibleRegion() );

    imIter.GoToBegin();

    long n = dataSz[0]*dataSz[1];
    for (long i = 0; i < n; ++i)
      {
        array2->set(i, imIter.Get());

        ++imIter;
      }

    typename ImageType::SpacingType sp = image->GetSpacing();
    array2->setSpacingX(sp[0]);
    array2->setSpacingY(sp[1]);

    return array2;
  }

  /*============================================================
    itkImageToArray3 */
  template< typename T >
  typename cArray3D< T >::PointerType itkImageToArray3(typename itk::Image< T, 3 >::Pointer image)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D
    typename ImageType::RegionType::SizeType dataSz = (image->GetLargestPossibleRegion()).GetSize();

    typename cArray3D<T>::PointerType array3( new cArray3D< T >(dataSz[0], dataSz[1], dataSz[2]) );

    typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;
    ImageIteratorType imIter(image, image->GetLargestPossibleRegion() );

    imIter.GoToBegin();

    long n = dataSz[0]*dataSz[1]*dataSz[2];
    for (long i = 0; i < n; ++i)
      {
        array3->set(i, imIter.Get());

        ++imIter;
      }

    typename ImageType::SpacingType sp = image->GetSpacing();
    array3->setSpacingX(sp[0]);
    array3->setSpacingY(sp[1]);
    array3->setSpacingZ(sp[2]);

    return array3;
  }


  /*============================================================
    cArray2ToItkImage  */
  template< typename T >
  typename itk::Image< T, 2 >::Pointer cArray2ToItkImage(typename cArray2D< T >::PointerType array2)
  {
    typedef itk::Image< T, 2 > ImageType; // 2 means 2D

    typename ImageType::Pointer image = ImageType::New();

    typename ImageType::IndexType imStart;
    imStart[0] = 0;  // first index on X
    imStart[1] = 0;  // first index on Y

    typename ImageType::SizeType  imSize;
    imSize[0] = array2->getSizeX();  // size along X
    imSize[1] = array2->getSizeY();  // size along Y

    typename ImageType::RegionType imRegion;
    imRegion.SetSize( imSize );
    imRegion.SetIndex( imStart );


    typename ImageType::SpacingType sp;
    sp[0] = array2->getSpacingX();
    sp[1] = array2->getSpacingY();

    image->SetSpacing(sp);
    image->SetRegions( imRegion );
    image->Allocate();

    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    IteratorType imIt(image, image->GetLargestPossibleRegion() );
    imIt.GoToBegin();

    for (long i = 0; !imIt.IsAtEnd(); ++imIt)
      {
        imIt.Set((*array2)[i]);
        ++i;
      }


    return image;
  }


  /*============================================================
    cArray3ToItkImage  */
  template< typename T >
  typename itk::Image< T, 3 >::Pointer cArray3ToItkImage(typename cArray3D< T >::PointerType array3)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D

    typename ImageType::Pointer image = ImageType::New();

    typename ImageType::IndexType imStart;
    imStart[0] = 0;  // first index on X
    imStart[1] = 0;  // first index on Y
    imStart[2] = 0;  // first index on Y

    typename ImageType::SizeType  imSize;
    imSize[0] = array3->getSizeX();  // size along X
    imSize[1] = array3->getSizeY();  // size along Y
    imSize[2] = array3->getSizeZ();  // size along Y

    typename ImageType::RegionType imRegion;
    imRegion.SetSize( imSize );
    imRegion.SetIndex( imStart );

    typename ImageType::SpacingType sp;
    sp[0] = array3->getSpacingX();
    sp[1] = array3->getSpacingY();
    sp[2] = array3->getSpacingZ();

    image->SetSpacing(sp);
    image->SetRegions( imRegion );
    image->Allocate();

    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    IteratorType imIt(image, image->GetLargestPossibleRegion() );
    imIt.GoToBegin();

    for (long i = 0; !imIt.IsAtEnd(); ++imIt)
      {
        imIt.Set((*array3)[i]);
        ++i;
      }


    return image;
  }



  /*============================================================
    cArray3ToItkImage  */
  template< typename T >
  typename itk::Image< T, 3 >::Pointer cArray3ToItkImage(const cArray3D< T >* array3)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D

    typename ImageType::Pointer image = ImageType::New();

    typename ImageType::IndexType imStart;
    imStart[0] = 0;  // first index on X
    imStart[1] = 0;  // first index on Y
    imStart[2] = 0;  // first index on Y

    typename ImageType::SizeType  imSize;
    imSize[0] = array3->getSizeX();  // size along X
    imSize[1] = array3->getSizeY();  // size along Y
    imSize[2] = array3->getSizeZ();  // size along Y

    typename ImageType::RegionType imRegion;
    imRegion.SetSize( imSize );
    imRegion.SetIndex( imStart );

    typename ImageType::SpacingType sp;
    sp[0] = array3->getSpacingX();
    sp[1] = array3->getSpacingY();
    sp[2] = array3->getSpacingZ();

    image->SetSpacing(sp);
    image->SetRegions( imRegion );
    image->Allocate();

    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    IteratorType imIt(image, image->GetLargestPossibleRegion() );
    imIt.GoToBegin();

    for (long i = 0; !imIt.IsAtEnd(); ++imIt)
      {
        imIt.Set((*array3)[i]);
        ++i;
      }


    return image;
  }

  /*============================================================
    crop array  */
  template< typename T >
  typename cArray2D< T >::PointerType cropArray(typename cArray2D< T >::PointerType array2, \
                                                long xMin, long xMax, long yMin, long yMax)
  {
    /*
      xMin is the min index, in the original array, to be *included* in
      the cropped array. While xMax is the max, in the original array,
      to be *included* in the cropped array.

      So the limit of xMin is 0, indicating: no crop at head. And the
      limit of xMax is nx-1, indicating: no crop at end.

      Same for y-dim.
    */
    long nx = array2->getSizeX();
    long ny = array2->getSizeY();

    if (xMin < 0 || xMax > nx-1 || yMin < 0 || yMax > ny-1 )
      {
        std::cerr<<"cropping out of bounds, exit.\n";
        raise(SIGABRT);
      }

    if (xMin > xMax || yMin > yMax)
      {
        std::cerr<<"xMin > xMax || yMin > yMax, exit.\n";
        raise(SIGABRT);
      }

    long szX = xMax - xMin + 1;
    long szY = yMax - yMin + 1;

    cArray2D< T >* cropped = new cArray2D< T >(szX, szY);

    for (long ix = xMin; ix <= xMax; ++ix)
      {
        for (long iy = yMin; iy <= yMax; ++iy)
          {
            // Notice the "=" in the upper limit, since we want to
            // include the xMax and yMax
            long iix = ix - xMin;
            long iiy = iy - yMin;

            (*cropped)(iix, iiy) = (*array2)(ix, iy);
          }
      }


    return cropped;
  }



  /* ============================================================
     cropArray 3D */
  template< typename T >
  typename cArray3D< T >::PointerType cropArray(typename cArray3D< T >::PointerType array3, \
                                                long xMin, long xMax, long yMin, long yMax, long zMin, long zMax)
  {
    /*
      xMin is the min index, in the original array, to be *included* in
      the cropped array. While xMax is the max, in the original array,
      to be *included* in the cropped array.

      So the limit of xMin is 0, indicating: no crop at head. And the
      limit of xMax is nx-1, indicating: no crop at end.

      Same for y-dim.
    */
    long nx = array3->getSizeX();
    long ny = array3->getSizeY();
    long nz = array3->getSizeZ();

    if (xMin < 0 || xMax > nx-1 || yMin < 0 || yMax > ny-1 || zMin < 0 || zMax > nz-1)
      {
        std::cerr<<"cropping out of bounds, exit.\n";
        raise(SIGABRT);
      }

    if (xMin > xMax || yMin > yMax || zMin > zMax)
      {
        std::cerr<<"xMin > xMax || yMin > yMax || zMin > zMax, exit.\n";
        raise(SIGABRT);
      }

    long szX = xMax - xMin + 1;
    long szY = yMax - yMin + 1;
    long szZ = zMax - zMin + 1;

    cArray3D< T >* cropped = new cArray3D< T >(szX, szY, szZ);

    for (long ix = xMin; ix <= xMax; ++ix)
      {
        for (long iy = yMin; iy <= yMax; ++iy)
          {
            for (long iz = zMin; iz <= zMax; ++iz)
              {
                // Notice the "=" in the upper limit, since we want to
                // include the xMax and yMax and zMax
                long iix = ix - xMin;
                long iiy = iy - yMin;
                long iiz = iz - zMin;

                (*cropped)(iix, iiy, iiz) = (*array3)(ix, iy, iz);
              }
          }
      }


    return cropped;
  }



  /*============================================================
    const pad array 2D */
  template< typename T >
  typename cArray2D< T >::PointerType constPadArray(typename cArray2D< T >::PointerType array2, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum)
  {
    long nx = array2->getSizeX();
    long ny = array2->getSizeY();

    long newNx = nx + xHeadNum + xTailNum;
    long newNy = ny + yHeadNum + yTailNum;

    typename cArray2D< T >::PointerType padded( new cArray2D< T >(newNx, newNy, padValue) );

    for (long ix = 0; ix < nx; ++ix)
      {
        for (long iy = 0; iy < ny; ++iy)
          {
            long newIx = ix + xHeadNum;
            long newIy = iy + yHeadNum;

            (*padded)(newIx, newIy) = (*array2)(ix, iy);
          }
      }

    return padded;
  }

  /*============================================================
    const pad array 3D  */
  template< typename T >
  typename cArray3D< T >::PointerType constPadArray(typename cArray3D< T >::PointerType array3, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum, \
                                                    long zHeadNum, long zTailNum)
  {
    long nx = array3->getSizeX();
    long ny = array3->getSizeY();
    long nz = array3->getSizeZ();

    long newNx = nx + xHeadNum + xTailNum;
    long newNy = ny + yHeadNum + yTailNum;
    long newNz = nz + zHeadNum + zTailNum;

    typename cArray3D< T >::PointerType padded( new cArray3D< T >(newNx, newNy, newNz, padValue) );

    for (long ix = 0; ix < nx; ++ix)
      {
        for (long iy = 0; iy < ny; ++iy)
          {
            for (long iz = 0; iz < nz; ++iz)
              {
                long newIx = ix + xHeadNum;
                long newIy = iy + yHeadNum;
                long newIz = iz + zHeadNum;

                (*padded)(newIx, newIy, newIz) = (*array3)(ix, iy, iz);
              }
          }
      }

    return padded;
  }



} // gth818n

////////////////////////////////////////////////////////////
#endif
