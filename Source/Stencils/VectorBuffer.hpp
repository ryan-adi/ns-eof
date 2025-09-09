#pragma once

#include "Definitions.hpp"

namespace Stencils {

  class VectorBuffer {
  public:
    RealType* u;
    RealType* v;
    RealType* w;
    VectorBuffer(int size, int dim);
    ~VectorBuffer();
  };
} // namespace Stencils