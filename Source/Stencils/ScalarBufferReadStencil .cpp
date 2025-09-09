#include "StdAfx.hpp"

#include "ScalarBufferReadStencil.hpp"
Stencils::ScalarBufferReadStencil::ScalarBufferReadStencil(const Parameters& parameters, ScalarField& scalarField):
  BufferStencil(parameters),
  scalarField_(scalarField),
  leftBuffer_(parameters.geometry.dim == 3 ? new RealType[2 * Ny_ * Nz_] : new RealType[2 * Ny_]),
  rightBuffer_(parameters.geometry.dim == 3 ? new RealType[Ny_ * Nz_] : new RealType[Ny_]),
  bottomBuffer_(parameters.geometry.dim == 3 ? new RealType[2 * Nx_ * Nz_] : new RealType[2 * Nx_]),
  topBuffer_(parameters.geometry.dim == 3 ? new RealType[Nx_ * Nz_] : new RealType[Nx_]),
  frontBuffer_(parameters.geometry.dim == 3 ? new RealType[2 * Nx_ * Ny_] : nullptr),
  backBuffer_(parameters.geometry.dim == 3 ? new RealType[Nx_ * Ny_] : nullptr) {}

Stencils::ScalarBufferReadStencil::~ScalarBufferReadStencil() {
  delete[] leftBuffer_;
  delete[] rightBuffer_;
  delete[] bottomBuffer_;
  delete[] topBuffer_;
  delete[] frontBuffer_;
  delete[] backBuffer_;
}

void Stencils::ScalarBufferReadStencil::applyLeftWall(int i, int j) {
  scalarField_.getScalar(i + 1, j) = leftBuffer_[j];
  scalarField_.getScalar(i, j)     = leftBuffer_[j + Ny_];
}
void Stencils::ScalarBufferReadStencil::applyRightWall(int i, int j) { scalarField_.getScalar(i, j) = rightBuffer_[j]; }
void Stencils::ScalarBufferReadStencil::applyBottomWall(int i, int j) {
  scalarField_.getScalar(i, j + 1) = bottomBuffer_[i];
  scalarField_.getScalar(i, j)     = bottomBuffer_[i + Nx_];
}
void Stencils::ScalarBufferReadStencil::applyTopWall(int i, int j) { scalarField_.getScalar(i, j) = topBuffer_[i]; }

void Stencils::ScalarBufferReadStencil::applyLeftWall(int i, int j, int k) {
  scalarField_.getScalar(i + 1, j, k) = leftBuffer_[j * Nz_ + k];
  scalarField_.getScalar(i, j, k)     = leftBuffer_[j * Nz_ + k + Ny_ * Nz_];
}
void Stencils::ScalarBufferReadStencil::applyRightWall(int i, int j, int k) { scalarField_.getScalar(i, j, k) = rightBuffer_[j * Nz_ + k]; }
void Stencils::ScalarBufferReadStencil::applyBottomWall(int i, int j, int k) {
  scalarField_.getScalar(i, j + 1, k) = bottomBuffer_[i * Nz_ + k];
  scalarField_.getScalar(i, j, k)     = bottomBuffer_[i * Nz_ + k + Nx_ * Nz_];
}
void Stencils::ScalarBufferReadStencil::applyTopWall(int i, int j, int k) { scalarField_.getScalar(i, j, k) = topBuffer_[i * Nz_ + k]; }
void Stencils::ScalarBufferReadStencil::applyFrontWall(int i, int j, int k) {
  scalarField_.getScalar(i, j, k + 1) = frontBuffer_[i * Ny_ + j];
  scalarField_.getScalar(i, j, k)     = frontBuffer_[i * Ny_ + j + Nx_ * Ny_];
}
void Stencils::ScalarBufferReadStencil::applyBackWall(int i, int j, int k) { scalarField_.getScalar(i, j, k) = backBuffer_[i * Ny_ + j]; }

RealType* Stencils::ScalarBufferReadStencil::getLeftBuffer() { return leftBuffer_; }
RealType* Stencils::ScalarBufferReadStencil::getRightBuffer() { return rightBuffer_; }
RealType* Stencils::ScalarBufferReadStencil::getBottomBuffer() { return bottomBuffer_; }
RealType* Stencils::ScalarBufferReadStencil::getTopBuffer() { return topBuffer_; }
RealType* Stencils::ScalarBufferReadStencil::getFrontBuffer() { return frontBuffer_; }
RealType* Stencils::ScalarBufferReadStencil::getBackBuffer() { return backBuffer_; }
