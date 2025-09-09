#include "StdAfx.hpp"

#include "ScalarBufferFillStencil.hpp"

Stencils::ScalarBufferFillStencil::ScalarBufferFillStencil(const Parameters& parameters, ScalarField& scalarField):
  BufferStencil(parameters),
  scalarField_(scalarField),
  leftBuffer_(parameters.geometry.dim == 3 ? new RealType[Ny_ * Nz_] : new RealType[Ny_]),
  rightBuffer_(parameters.geometry.dim == 3 ? new RealType[2 * Ny_ * Nz_] : new RealType[2 * Ny_]),
  bottomBuffer_(parameters.geometry.dim == 3 ? new RealType[Nx_ * Nz_] : new RealType[Nx_]),
  topBuffer_(parameters.geometry.dim == 3 ? new RealType[2 * Nx_ * Nz_] : new RealType[2 * Nx_]),
  frontBuffer_(parameters.geometry.dim == 3 ? new RealType[Nx_ * Ny_] : nullptr),
  backBuffer_(parameters.geometry.dim == 3 ? new RealType[2 * Nx_ * Ny_] : nullptr) {}

Stencils::ScalarBufferFillStencil::~ScalarBufferFillStencil() {
  delete[] leftBuffer_;
  delete[] rightBuffer_;
  delete[] bottomBuffer_;
  delete[] topBuffer_;
  delete[] frontBuffer_;
  delete[] backBuffer_;
}

void Stencils::ScalarBufferFillStencil::applyLeftWall(int i, int j) { leftBuffer_[j] = scalarField_.getScalar(i + 2, j); }
void Stencils::ScalarBufferFillStencil::applyRightWall(int i, int j) {
  rightBuffer_[j]       = scalarField_.getScalar(i - 1, j);
  rightBuffer_[j + Ny_] = scalarField_.getScalar(i - 2, j);
}
void Stencils::ScalarBufferFillStencil::applyBottomWall(int i, int j) { bottomBuffer_[i] = scalarField_.getScalar(i, j + 2); }
void Stencils::ScalarBufferFillStencil::applyTopWall(int i, int j) {
  topBuffer_[i]       = scalarField_.getScalar(i, j - 1);
  topBuffer_[i + Nx_] = scalarField_.getScalar(i, j - 2);
}
void Stencils::ScalarBufferFillStencil::applyLeftWall(int i, int j, int k) { leftBuffer_[j * Nz_ + k] = scalarField_.getScalar(i + 2, j, k); }
void Stencils::ScalarBufferFillStencil::applyRightWall(int i, int j, int k) {
  rightBuffer_[j * Nz_ + k]             = scalarField_.getScalar(i - 1, j, k);
  rightBuffer_[j * Nz_ + k + Ny_ * Nz_] = scalarField_.getScalar(i - 2, j, k);
}
void Stencils::ScalarBufferFillStencil::applyBottomWall(int i, int j, int k) { bottomBuffer_[i * Nz_ + k] = scalarField_.getScalar(i, j + 2, k); }
void Stencils::ScalarBufferFillStencil::applyTopWall(int i, int j, int k) {
  topBuffer_[i * Nz_ + k]             = scalarField_.getScalar(i, j - 1, k);
  topBuffer_[i * Nz_ + k + Nx_ * Nz_] = scalarField_.getScalar(i, j - 2, k);
}
void Stencils::ScalarBufferFillStencil::applyFrontWall(int i, int j, int k) { frontBuffer_[i * Ny_ + j] = scalarField_.getScalar(i, j, k + 2); }
void Stencils::ScalarBufferFillStencil::applyBackWall(int i, int j, int k) {
  backBuffer_[i * Ny_ + j]             = scalarField_.getScalar(i, j, k - 1);
  backBuffer_[i * Ny_ + j + Nx_ * Ny_] = scalarField_.getScalar(i, j, k - 2);
}

RealType* Stencils::ScalarBufferFillStencil::getLeftBuffer() { return leftBuffer_; }
RealType* Stencils::ScalarBufferFillStencil::getRightBuffer() { return rightBuffer_; }
RealType* Stencils::ScalarBufferFillStencil::getBottomBuffer() { return bottomBuffer_; }
RealType* Stencils::ScalarBufferFillStencil::getTopBuffer() { return topBuffer_; }
RealType* Stencils::ScalarBufferFillStencil::getFrontBuffer() { return frontBuffer_; }
RealType* Stencils::ScalarBufferFillStencil::getBackBuffer() { return backBuffer_; }
