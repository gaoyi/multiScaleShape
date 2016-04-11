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

#ifndef cArray2D_txx_
#define cArray2D_txx_

#include "cArray2D.h"

#include <iostream>

namespace gth818n
{
  template<typename T>
  cArray2D<T>::cArray2D()
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;

    m_data = 0;

    m_nx = 0;
    m_ny = 0;
  }


  template<typename T>
  cArray2D<T>::~cArray2D()
  {
    clear();
  }


  template<typename T>
  cArray2D<T>::cArray2D(long nx, long ny)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;

    m_nx = nx;
    m_ny = ny;

    m_data = new T[nx*ny];
  }

  template<typename T>
  cArray2D<T>::cArray2D(long nx, long ny, T initVal)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;

    m_nx = nx;
    m_ny = ny;

    m_data = new T[nx*ny];

    fill(initVal);
  }


  template<typename T>
  cArray2D<T>::cArray2D(const cArray2D& a)
  {
    m_spacingX = a.m_spacingX;
    m_spacingY = a.m_spacingY;

    m_nx = a.m_nx;
    m_ny = a.m_ny;

    long n = m_nx*m_ny;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = a.m_data[i];
      }
  }

  template<typename T>
  cArray2D<T>::cArray2D(const PointerType p)
  {
    m_spacingX = p->m_spacingX;
    m_spacingY = p->m_spacingY;

    m_nx = p->m_nx;
    m_ny = p->m_ny;

    long n = m_nx*m_ny;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = p->m_data[i];
      }
  }

  template<typename T>
  cArray2D< T >& cArray2D< T >::operator=(cArray2D< T > const& rhs)
  {
    if (this == &rhs)
      {
        return *this;
      }

    // this != &rhs
    if (rhs.m_data)
      {
        // copy nx, ny
        if (this->m_nx != rhs.m_nx || this->m_ny != rhs.m_ny)
          {
            this->setSize(rhs.m_nx, rhs.m_ny);
          }

        // copy data
        long n = m_nx*m_ny;
        for (long i = 0; i <= n-1; i++)
          {
            this->m_data[i] = rhs.m_data[i];
          }

        // copy spacing
        m_spacingX = rhs.m_spacingX;
        m_spacingY = rhs.m_spacingY;
      }

    return *this;
  }



  template<typename T>
  void cArray2D<T>::clear()
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

        m_nx = 0;
        m_ny = 0;

        m_spacingX = 0;
        m_spacingY = 0;
      }

    return;
  }


  template<typename T>
  bool cArray2D<T>::setSize(long nx, long ny)
  {
     if (this->m_data)
       {
         // m_data already allocated
         if ( (this->m_nx == nx) && (this->m_ny == ny) )
           {
             // if no change in size, do not reallocate.
             return false;
           }

         this->clear();

         this->m_nx = nx;
         this->m_ny = ny;

         this->m_data = new T[nx*ny];
       }
     else
       {
         this->m_nx = nx;
         this->m_ny = ny;

         this->m_data = new T[nx*ny];
       }

     return true;
  }


  template<typename T>
  void cArray2D<T>::fill(T val)
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }
    else
      {
        long n = m_nx*m_ny;
        for (long i = 0; i <= n-1; ++i)
          {
            m_data[i] = val;
          }
      }

    return;
  }

} //gth818n

#endif
