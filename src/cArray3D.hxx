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


#ifndef cArray3D_txx_
#define cArray3D_txx_

#include "cArray3D.h"

#include <iostream>

namespace gth818n
{
  template<typename T>
  cArray3D<T>::cArray3D()
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;

    m_data = 0;

    m_nx = 0;
    m_ny = 0;
    m_nz = 0;

    m_strideX = 1;
    m_strideY = 1;
    m_strideZ = 1;
  }


  template<typename T>
  cArray3D<T>::~cArray3D()
  {
    clear();
  }


  template<typename T>
  cArray3D<T>::cArray3D(long nx, long ny, long nz)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;

    m_nx = nx;
    m_ny = ny;
    m_nz = nz;

    m_strideX = 1;
    m_strideY = m_nx;
    m_strideZ = m_nx*m_ny;

    m_data = new T[nx*ny*nz];
  }

  template<typename T>
  cArray3D<T>::cArray3D(long nx, long ny, long nz, T initVal)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;

    m_nx = nx;
    m_ny = ny;
    m_nz = nz;

    m_strideX = 1;
    m_strideY = m_nx;
    m_strideZ = m_nx*m_ny;

    m_data = new T[nx*ny*nz];

    fill(initVal);
  }


  template<typename T>
  cArray3D<T>::cArray3D(const cArray3D &a)
  {
    m_spacingX = a.m_spacingX;
    m_spacingY = a.m_spacingY;
    m_spacingZ = a.m_spacingZ;

    m_nx = a.m_nx;
    m_ny = a.m_ny;
    m_nz = a.m_nz;

    long n = m_nx*m_ny*m_nz;

    m_data = new T[n];

    m_strideX = 1;
    m_strideY = m_nx;
    m_strideZ = m_nx*m_ny;


    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = a.m_data[i];
      }
  }



  template<typename T>
  cArray3D<T>::cArray3D(const PointerType p)
  {
    m_spacingX = p->m_spacingX;
    m_spacingY = p->m_spacingY;
    m_spacingZ = p->m_spacingZ;

    m_nx = p->m_nx;
    m_ny = p->m_ny;
    m_nz = p->m_nz;

    m_strideX = 1;
    m_strideY = m_nx;
    m_strideZ = m_nx*m_ny;


    long n = m_nx*m_ny*m_nz;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = p->m_data[i];
      }
  }


  template<typename T>
  cArray3D< T > &cArray3D< T >::operator=(cArray3D< T > const &rhs)
  {
    if (this == &rhs)
      {
        return *this;
      }

    // this != &rhs
    if (rhs.m_data)
      {
        // copy nx, ny, nz
        if (this->m_nx != rhs.m_nx || this->m_ny != rhs.m_ny || this->m_nz != rhs.m_nz)
          {
            this->setSize(rhs.m_nx, rhs.m_ny, rhs.m_nz);
          }

        m_strideX = 1;
        m_strideY = m_nx;
        m_strideZ = m_nx*m_ny;

        // copy data
        long n = m_nx*m_ny*m_nz;
        for (long i = 0; i <= n-1; i++)
          {
            this->m_data[i] = rhs.m_data[i];
          }

        // copy spacing
        m_spacingX = rhs.m_spacingX;
        m_spacingY = rhs.m_spacingY;
        m_spacingZ = rhs.m_spacingZ;
      }

    return *this;
  }


  template<typename T>
  void cArray3D<T>::clear()
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }
    else
      {
        delete[] m_data;
        m_data = 0;
      }

    return;
  }

  template<typename T>
  bool cArray3D<T>::setSize(long nx, long ny, long nz)
  {
     if (this->m_data)
       {
         // m_data already allocated
         if ( (this->m_nx == nx) && (this->m_ny == ny) && (this->m_nz == nz) )
           {
             // if no change in size, do not reallocate.
             return false;
           }

         this->clear();

         this->m_nx = nx;
         this->m_ny = ny;
         this->m_nz = nz;

         this->m_data = new T[nx*ny*nz];
       }
     else
       {
         this->m_nx = nx;
         this->m_ny = ny;
         this->m_nz = nz;

         this->m_data = new T[nx*ny*nz];
       }

     return true;
  }


  template<typename T>
  void cArray3D<T>::fill(T val)
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }
    else
      {
        long n = m_nx*m_ny*m_nz;
        for (long i = 0; i <= n-1; ++i)
          {
            m_data[i] = val;
          }
      }

    return;
  }

  template<typename T>
  T cArray3D< T >::max_value() const
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }

    T maxValue = m_data[0];
    long n = m_nx*m_ny*m_nz;
    for (long i = 0; i <= n-1; ++i)
      {
        T currentValue = m_data[i];
        maxValue = maxValue>currentValue?maxValue:currentValue;
      }

    return maxValue;
  }

  template<typename T>
  T cArray3D< T >::min_value() const
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }

    T minValue = m_data[0];
    long n = m_nx*m_ny*m_nz;
    for (long i = 0; i <= n-1; ++i)
      {
        T currentValue = m_data[i];
        minValue = minValue>currentValue?currentValue:minValue;
      }

    return minValue;
  }

  template<typename T>
  void cArray3D< T >::extreme_values(T &minValue, T &maxValue) const
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }

    minValue = m_data[0];
    maxValue = minValue;

    long n = m_nx*m_ny*m_nz;
    for (long i = 0; i <= n-1; ++i)
      {
        T currentValue = m_data[i];
        maxValue = maxValue>currentValue?maxValue:currentValue;
        minValue = minValue>currentValue?currentValue:minValue;
      }

    return;
  }



} //namespace

#endif
