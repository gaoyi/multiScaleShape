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


/*********************************************************************************
 * Open nrrd file
 * dwt
 * wimage to band
 * band to wimage
 * idwt
 **********************************************************************************/

#include "wavelet.h"

#include "cArrayOp.h"


#include <ctime>

int main(int argc, char* argv[])
{
  if (argc < 3)
    {
      std::cerr<<"Usage: "<<argv[0]<<" inputImage outputImage\n";
      exit(-1);
    }

  const char* inputFileName = argv[1];
  const char* outputFileName = argv[2];

  typedef double pixel_t;
  typedef itk::Image<pixel_t, 3> itkImage_t;

  itkImage_t::Pointer im = gth818n::readImage<pixel_t, 3>(inputFileName);


  ////////////////////////////////////////////////////////////
  // begin timing
  double t1, t2;
  t1 = clock();

  itkImage_t::Pointer dwtCoef = gth818n::dwt3< double >(im);

  t2 = clock();
  t1 = (t2 - t1)/CLOCKS_PER_SEC;
  std::cerr<<"Ellapsed time of DWT = "<<t1<<std::endl;
  // end timing
  ////////////////////////////////////////////////////////////


  /**
   * Wimage to band
   */
  boost::shared_ptr< std::vector< vnl_vector< double > > > bands = gth818n::wImageToBandStructure<double>(dwtCoef);

  /**
   * print band info
   */
  std::cout<<"num of bands = "<<bands->size()<<std::endl;

  for (unsigned long ib = 0; ib < bands->size(); ++ib )
    {
      const vnl_vector< double >& thisBand = (*bands)[ib];
      std::cout<<"In band "<<ib<<", has "<<thisBand.size()<<" coefs.\n";
    }




  /**
   * band to Wimage
   */
  long nx = im->GetLargestPossibleRegion().GetSize()[0];
  long ny = im->GetLargestPossibleRegion().GetSize()[1];
  long nz = im->GetLargestPossibleRegion().GetSize()[2];
  itkImage_t::Pointer newWimage = gth818n::bandStructureToWimage<double>( bands, nx, ny, nz);

  // inverse WT
  ////////////////////////////////////////////////////////////
  // begin timing
  t1 = clock();

  itkImage_t::Pointer reconIm = gth818n::idwt3< double >(newWimage);

  t2 = clock();
  t1 = (t2 - t1)/CLOCKS_PER_SEC;
  std::cerr<<"Ellapsed time of IDWT = "<<t1<<std::endl;
  // end timing
  ////////////////////////////////////////////////////////////

  gth818n::writeImage< pixel_t, 3 >(reconIm, outputFileName);

  return 0;
}
