#pragma once

#include "Parameters.hpp"

namespace Stencils {

  /** Interface for operations on the global (or parallel) boundary
   * New definition of class, equivalent to BoundaryStencil, better suited for our parallel algorithm
   */
  class BufferStencil {
  protected:
    const Parameters& parameters_;
    const int         Nx_, Ny_, Nz_;

  public:
    BufferStencil(const Parameters& parameters):
      parameters_(parameters),
      Nx_(parameters.parallel.localSize[0] + 3),
      Ny_(parameters.parallel.localSize[1] + 3),
      Nz_(parameters.parallel.localSize[2] + 3) {}

    virtual ~BufferStencil() = default;

    virtual void applyLeftWall(int i, int j)          = 0;
    virtual void applyRightWall(int i, int j)         = 0;
    virtual void applyBottomWall(int i, int j)        = 0;
    virtual void applyTopWall(int i, int j)           = 0;
    virtual void applyLeftWall(int i, int j, int k)   = 0;
    virtual void applyRightWall(int i, int j, int k)  = 0;
    virtual void applyBottomWall(int i, int j, int k) = 0;
    virtual void applyTopWall(int i, int j, int k)    = 0;
    virtual void applyFrontWall(int i, int j, int k)  = 0;
    virtual void applyBackWall(int i, int j, int k)   = 0;

    virtual RealType* getLeftBuffer()   = 0;
    virtual RealType* getRightBuffer()  = 0;
    virtual RealType* getBottomBuffer() = 0;
    virtual RealType* getTopBuffer()    = 0;
    virtual RealType* getFrontBuffer()  = 0;
    virtual RealType* getBackBuffer()   = 0;
  };


} // namespace Stencils
