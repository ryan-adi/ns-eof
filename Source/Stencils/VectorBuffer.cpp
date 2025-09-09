#include "StdAfx.hpp"

#include "VectorBuffer.hpp"

Stencils::VectorBuffer::VectorBuffer(int size, int dim) {
  u = nullptr;
  v = nullptr;
  w = nullptr;
  if (size > 0) {
    if (dim == 2) {
      u = new RealType[2 * size];
      v = u + size;
    } else if (dim == 3) {
      u = new RealType[3 * size];
      v = u + size;
      w = v + size;
    }
  }
}

Stencils::VectorBuffer::~VectorBuffer() { delete[] u; }