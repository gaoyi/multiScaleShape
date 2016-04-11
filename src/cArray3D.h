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


#ifndef cArray3D_h_
#define cArray3D_h_

#include "basicHeaders.h"

namespace gth818n
{
  template<typename T>
  class cArray3D
  {
  public:
    typedef cArray3D< T > Self;
    typedef boost::shared_ptr< Self > Pointer;
    typedef Pointer PointerType; // just for back compatability

    typedef T data_t;


    /* functions */
    cArray3D();
    cArray3D(long nx, long ny, long nz);
    cArray3D(long nx, long ny, long nz, T initVal);
    cArray3D(const cArray3D &a);
    cArray3D(const PointerType p);

    ~cArray3D();

    cArray3D< T > &operator=(cArray3D< T > const &rhs);


    inline T const & operator()(long x, long y, long z) const
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );

      return m_data[z*m_ny*m_nx + y*m_nx + x];
    }

    inline T& operator()(long x, long y, long z)
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );

      return m_data[z*m_ny*m_nx + y*m_nx + x];
    }


    inline T const & operator[](long idx) const
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz-1);

      return m_data[idx];
    }

    inline T& operator[](long idx)
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz-1);

      return m_data[idx];
    }

    //////////////////////////////////////////////////
    inline T get(long x, long y, long z) const
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );

      return m_data[z*m_ny*m_nx + y*m_nx + x];
    }

    inline T get(long idx) const
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz-1);

      return m_data[idx];
    }


    void set(long x, long y, long z, T v)
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );

       m_data[z*m_ny*m_nx + y*m_nx + x] = v;
    }


    void set(long idx, T v)
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz-1);

      m_data[idx] = v;
    }


    long getSizeX() const {return m_nx;}
    long getSizeY() const {return m_ny;}
    long getSizeZ() const {return m_nz;}

    inline long sub2idx(long ix, long iy, long iz) {return ix + iy*m_nx + iz*m_nx*m_ny;}


    double getSpacingX() const {return m_spacingX;}
    double getSpacingY() const {return m_spacingY;}
    double getSpacingZ() const {return m_spacingZ;}

    void setSpacingX(double spX) {m_spacingX = spX;}
    void setSpacingY(double spY) {m_spacingY = spY;}
    void setSpacingZ(double spZ) {m_spacingZ = spZ;}


    //    const T* getDataPointer() const {return m_data;}
    T* getDataPointer() const {return m_data;}

    T max_value() const;
    T min_value() const;
    void extreme_values(T &minValue, T &maxValue) const;

    bool setSize(long nx, long ny, long nz);


    void fill(T val);
    void clear();

    long m_strideX;
    long m_strideY;
    long m_strideZ;

    /* data */
  protected:
    T *m_data;
    long m_nx, m_ny, m_nz;
    double m_spacingX, m_spacingY, m_spacingZ;
  };

} //gth818n

#include "cArray3D.hxx"


#endif
/* definition in .hxx*/
