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

#ifndef cArray2D_h_
#define cArray2D_h_

#include "basicHeaders.h"

namespace gth818n
{
  template<typename T>
  class cArray2D
  {
  public:
    typedef cArray2D< T > Self;
    typedef boost::shared_ptr< Self > Pointer;
    typedef Pointer PointerType ; // just for back compatability

    /* functions */
    cArray2D();
    cArray2D(long nx, long ny);
    cArray2D(long nx, long ny, T initVal);
    cArray2D(const cArray2D& a);
    cArray2D(const PointerType p);
    ~cArray2D();

    cArray2D< T >& operator=(cArray2D< T > const& rhs);


    inline T const& operator()(long x, long y) const
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );

      return m_data[y*m_nx + x];
    }


    inline T& operator()(long x, long y)
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );

      return m_data[y*m_nx + x];
    }

    inline T& operator[](long idx)
    {
      assert(idx >= 0 && idx <= m_nx*m_ny-1);
      return m_data[idx];
    }

    inline T const& operator[](long idx) const
    {
      assert(idx >= 0 && idx <= m_nx*m_ny-1);
      return m_data[idx];
    }


    inline T get(long x, long y) const
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );

      return m_data[y*m_nx + x];
    }

    inline T get(long idx) const
    {
      assert(idx >= 0 && idx <= m_nx*m_ny-1);
      return m_data[idx];
    }


    void set(long x, long y, T v)
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );

      m_data[y*m_nx + x] = v;
    }


    void set(long idx, T v)
    {
      assert(idx >= 0 && idx <= m_nx*m_ny-1);
      m_data[idx] = v;
    }


    long getSizeX() const {return m_nx;}
    long getSizeY() const {return m_ny;}

    double getSpacingX() const {return m_spacingX;}
    double getSpacingY() const {return m_spacingY;}

    bool setSize(long nx, long ny);

    void setSpacingX(double spX) {m_spacingX = spX;}
    void setSpacingY(double spY) {m_spacingY = spY;}

    //    const T* getDataPointer() const {return m_data;}
    T* getDataPointer() const {return m_data;}

    void fill(T val);
    void clear();


  protected:
    /* data */
    T *m_data;
    long m_nx, m_ny;
    double m_spacingX, m_spacingY;
  };


} //gth818n

#include "cArray2D.hxx"


#endif
/* definition in .hxx*/
