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

#include "cArrayOp.h"


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

  typedef double pixel_t;
  typedef itk::Image<pixel_t, 3> itkImage_t;

  itkImage_t::Pointer im = gth818n::readImage<pixel_t, 3>(inputFileName);


  ////////////////////////////////////////////////////////////
  // begin timing
  double t1, t2;
  t1 = clock();

  itkImage_t::Pointer dwtCoef = gth818n::dwt3< double >(im);

  gth818n::writeImage< pixel_t, 3 >(dwtCoef, "coef.nrrd");


  t2 = clock();
  t1 = (t2 - t1)/CLOCKS_PER_SEC;
  std::cerr<<"Ellapsed time of DWT = "<<t1<<std::endl;
  // end timing
  ////////////////////////////////////////////////////////////




  // inverse WT
  ////////////////////////////////////////////////////////////
  // begin timing
  t1 = clock();

  itkImage_t::Pointer reconIm = gth818n::idwt3< double >(dwtCoef);

  t2 = clock();
  t1 = (t2 - t1)/CLOCKS_PER_SEC;
  std::cerr<<"Ellapsed time of IDWT = "<<t1<<std::endl;
  // end timing
  ////////////////////////////////////////////////////////////


  gth818n::writeImage< pixel_t, 3 >(reconIm, outputFileName);





  return 0;
}
