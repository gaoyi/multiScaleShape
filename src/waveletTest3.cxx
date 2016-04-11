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


#include "wavelet.h"

#include <ctime>

int main(int argc, char* argv[])
{
  if (argc != 4)
    {
      std::cerr<<"Usage: "<<argv[0]<<" inputImage outputImage thld\n";
      exit(-1);
    }

  const char* inputFileName = argv[1];
  const char* outputFileName = argv[2];
  double thld = atof(argv[3]);

  typedef double PixelType;
  typedef gth818n::cArray3D< PixelType > ArrayType;

  ArrayType::PointerType im( gth818n::readImageToArray3< PixelType >(inputFileName) );

  ArrayType::PointerType imRightSize = gth818n::padToRightSize<double>(im);

  im.reset();// free memory of im by dropping its ref count by 1


  ////////////////////////////////////////////////////////////
  // begin timing
  double t1, t2;
  t1 = clock();

  ArrayType::PointerType dwtCoef = gth818n::dwt3< double >(imRightSize);

  gth818n::saveAsImage3< PixelType >(dwtCoef, "coef.nrrd");

  long nx = dwtCoef->getSizeX();
  long ny = dwtCoef->getSizeY();
  long nz = dwtCoef->getSizeZ();


  t2 = clock();
  t1 = (t2 - t1)/CLOCKS_PER_SEC;
  std::cerr<<"Ellapsed time of DWT = "<<t1<<std::endl;
  // end timing
  ////////////////////////////////////////////////////////////


  gth818n::hardThresholdByPercentage< double >(dwtCoef, thld);

  gth818n::saveAsImage3< PixelType >(dwtCoef, "coefPost.nrrd");
  //  std::cout<<"reduction ratio = "<<k/(double)(nx*ny*nz)<<std::endl;

  typedef vnl_vector< double > VnlDoubleVectorType;

  typedef std::vector< VnlDoubleVectorType > WaveletBandsType;
  typedef boost::shared_ptr< WaveletBandsType > WaveletBandsPointerType;

  WaveletBandsPointerType bands = gth818n::coefArrayToBandStructure< double >(dwtCoef);

  std::cout<<"totally "<<bands->size()<<" bands.\n";


  ArrayType::PointerType reconCoef = gth818n::bandStructureToCoefArray3D< double >( bands, nx, ny, nz);
  gth818n::saveAsImage3< PixelType >(reconCoef, "coefPostRecon.nrrd");

  // inverse WT
  ////////////////////////////////////////////////////////////
  // begin timing
  t1 = clock();

  ArrayType::PointerType reconIm = gth818n::idwt3< double >(reconCoef);

  t2 = clock();
  t1 = (t2 - t1)/CLOCKS_PER_SEC;
  std::cerr<<"Ellapsed time of IDWT = "<<t1<<std::endl;
  // end timing
  ////////////////////////////////////////////////////////////


  gth818n::saveAsImage3< PixelType >(reconIm, outputFileName);


  return 0;
}
