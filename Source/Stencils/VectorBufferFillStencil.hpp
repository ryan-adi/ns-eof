#pragma once

#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Stencils/BufferStencil.hpp"
#include "Stencils/VectorBuffer.hpp"

namespace Stencils {

  class VectorBufferFillStencil: public BufferStencil {
  private:
    VectorField& vectorField_;

    VectorBuffer leftBuffer_;
    VectorBuffer rightBuffer_;
    VectorBuffer bottomBuffer_;
    VectorBuffer topBuffer_;
    VectorBuffer frontBuffer_;
    VectorBuffer backBuffer_;

  public:
    VectorBufferFillStencil(const Parameters& parameters, VectorField& vectorField);
    ~VectorBufferFillStencil() override;

    void applyLeftWall(int i, int j) override;
    void applyRightWall(int i, int j) override;
    void applyBottomWall(int i, int j) override;
    void applyTopWall(int i, int j) override;

    void applyLeftWall(int i, int j, int k) override;
    void applyRightWall(int i, int j, int k) override;
    void applyBottomWall(int i, int j, int k) override;
    void applyTopWall(int i, int j, int k) override;
    void applyFrontWall(int i, int j, int k) override;
    void applyBackWall(int i, int j, int k) override;

    RealType* getLeftBuffer() override;
    RealType* getRightBuffer() override;
    RealType* getBottomBuffer() override;
    RealType* getTopBuffer() override;
    RealType* getFrontBuffer() override;
    RealType* getBackBuffer() override;
  };

} // namespace Stencils