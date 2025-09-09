#pragma once

#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Stencils/BufferStencil.hpp"
namespace Stencils {

  class ScalarBufferFillStencil: public BufferStencil {
  private:
      ScalarField& scalarField_;

    RealType* leftBuffer_;
    RealType* rightBuffer_;
    RealType* bottomBuffer_;
    RealType* topBuffer_;
    RealType* frontBuffer_;
    RealType* backBuffer_;

  public:
    ScalarBufferFillStencil(const Parameters& parameters, ScalarField& scalarField);
    ~ScalarBufferFillStencil() override;

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