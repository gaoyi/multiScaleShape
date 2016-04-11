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
  typedef gth818n::cArray2D< PixelType > arrayType;

  arrayType::PointerType im( gth818n::readImageToArray2< PixelType >(inputFileName) );

  arrayType::PointerType imRightSize = gth818n::padToRightSize<double>(im);

  arrayType::PointerType dwtCoef = gth818n::dwt2< double >(imRightSize);

  gth818n::saveAsImage2< PixelType >(dwtCoef, "coef.nrrd");

  gth818n::hardThresholdByPercentage< double >(dwtCoef, thld);


  gth818n::saveAsImage2< PixelType >(dwtCoef, "coefPost.nrrd");

  typedef vnl_vector< double > VnlDoubleVectorType;

  typedef std::vector< VnlDoubleVectorType > WaveletBandsType;
  typedef boost::shared_ptr< WaveletBandsType > WaveletBandsPointerType;

  WaveletBandsPointerType bands = gth818n::coefArrayToBandStructure< double >(dwtCoef);


  std::cout<<"totally "<<bands->size()<<" bands.\n";

  long nx = dwtCoef->getSizeX();
  long ny = dwtCoef->getSizeY();
  arrayType::PointerType reconCoef = gth818n::bandStructureToCoefArray2D< double >( bands, nx, ny);
  gth818n::saveAsImage2< PixelType >(reconCoef, "coefPostRecon.nrrd");
  //  std::cout<<"reduction ratio = "<<k/(double)(nx*ny)<<std::endl;


  // inverse WT
  arrayType::PointerType reconIm = gth818n::idwt2< double >(reconCoef);

  gth818n::saveAsImage2< PixelType >(reconIm, outputFileName);


  return 0;
}
