#include "StdAfx.hpp"

#include "VectorBufferReadStencil.hpp"


Stencils::VectorBufferReadStencil::VectorBufferReadStencil(const Parameters& parameters, VectorField& vectorField):
  BufferStencil(parameters),
  vectorField_(vectorField),
  leftBuffer_(2 * (parameters.geometry.dim == 3 ? Ny_ * Nz_ : Ny_), parameters.geometry.dim),
  rightBuffer_((parameters.geometry.dim == 3 ? Ny_ * Nz_ : Ny_), parameters.geometry.dim),
  bottomBuffer_(2 * (parameters.geometry.dim == 3 ? Nx_ * Nz_ : Nx_), parameters.geometry.dim),
  topBuffer_((parameters.geometry.dim == 3 ? Nx_ * Nz_ : Nx_), parameters.geometry.dim),
  frontBuffer_(2 * (parameters.geometry.dim == 3 ? Nx_ * Ny_ : 0), parameters.geometry.dim),
  backBuffer_((parameters.geometry.dim == 3 ? Nx_ * Ny_ : 0), parameters.geometry.dim) {}

Stencils::VectorBufferReadStencil::~VectorBufferReadStencil() {}

void Stencils::VectorBufferReadStencil::applyLeftWall(int i, int j) {
  RealType* velocity_inner = vectorField_.getVector(i + 1, j);
  velocity_inner[0]        = leftBuffer_.u[j];
  velocity_inner[1]        = leftBuffer_.v[j];
  RealType* velocity_outer = vectorField_.getVector(i, j);
  velocity_outer[0]        = leftBuffer_.u[j + Ny_];
  velocity_outer[1]        = leftBuffer_.v[j + Ny_];
}
void Stencils::VectorBufferReadStencil::applyRightWall(int i, int j) {
  RealType* velocity = vectorField_.getVector(i, j);
  velocity[0]        = rightBuffer_.u[j];
  velocity[1]        = rightBuffer_.v[j];
}
void Stencils::VectorBufferReadStencil::applyBottomWall(int i, int j) {
  RealType* velocity_inner = vectorField_.getVector(i, j + 1);
  velocity_inner[0]        = bottomBuffer_.u[i];
  velocity_inner[1]        = bottomBuffer_.v[i];
  RealType* velocity_outer = vectorField_.getVector(i, j);
  velocity_outer[0]        = bottomBuffer_.u[i + Nx_];
  velocity_outer[1]        = bottomBuffer_.v[i + Nx_];
}
void Stencils::VectorBufferReadStencil::applyTopWall(int i, int j) {
  RealType* velocity = vectorField_.getVector(i, j);
  velocity[0]        = topBuffer_.u[i];
  velocity[1]        = topBuffer_.v[i];
}

void Stencils::VectorBufferReadStencil::applyLeftWall(int i, int j, int k) {
  RealType* velocity_inner = vectorField_.getVector(i + 1, j, k);
  velocity_inner[0]        = leftBuffer_.u[j * Nz_ + k];
  velocity_inner[1]        = leftBuffer_.v[j * Nz_ + k];
  velocity_inner[2]        = leftBuffer_.w[j * Nz_ + k];
  RealType* velocity_outer = vectorField_.getVector(i, j, k);
  velocity_outer[0]        = leftBuffer_.u[j * Nz_ + k + Ny_ * Nz_];
  velocity_outer[1]        = leftBuffer_.v[j * Nz_ + k + Ny_ * Nz_];
  velocity_outer[2]        = leftBuffer_.w[j * Nz_ + k + Ny_ * Nz_];
}
void Stencils::VectorBufferReadStencil::applyRightWall(int i, int j, int k) {
  RealType* velocity = vectorField_.getVector(i, j, k);
  velocity[0]        = rightBuffer_.u[j * Nz_ + k];
  velocity[1]        = rightBuffer_.v[j * Nz_ + k];
  velocity[2]        = rightBuffer_.w[j * Nz_ + k];
}
void Stencils::VectorBufferReadStencil::applyBottomWall(int i, int j, int k) {
  RealType* velocity_inner = vectorField_.getVector(i, j + 1, k);
  velocity_inner[0]        = bottomBuffer_.u[i * Nz_ + k];
  velocity_inner[1]        = bottomBuffer_.v[i * Nz_ + k];
  velocity_inner[2]        = bottomBuffer_.w[i * Nz_ + k];
  RealType* velocity_outer = vectorField_.getVector(i, j, k);
  velocity_outer[0]        = bottomBuffer_.u[i * Nz_ + k + Nx_ * Nz_];
  velocity_outer[1]        = bottomBuffer_.v[i * Nz_ + k + Nx_ * Nz_];
  velocity_outer[2]        = bottomBuffer_.w[i * Nz_ + k + Nx_ * Nz_];
}
void Stencils::VectorBufferReadStencil::applyTopWall(int i, int j, int k) {
  RealType* velocity = vectorField_.getVector(i, j, k);
  velocity[0]        = topBuffer_.u[i * Nz_ + k];
  velocity[1]        = topBuffer_.v[i * Nz_ + k];
  velocity[2]        = topBuffer_.w[i * Nz_ + k];
}
void Stencils::VectorBufferReadStencil::applyFrontWall(int i, int j, int k) {
  RealType* velocity_inner = vectorField_.getVector(i, j, k + 1);
  velocity_inner[0]        = frontBuffer_.u[i * Ny_ + j];
  velocity_inner[1]        = frontBuffer_.v[i * Ny_ + j];
  velocity_inner[2]        = frontBuffer_.w[i * Ny_ + j];
  RealType* velocity_outer = vectorField_.getVector(i, j, k);
  velocity_outer[0]        = frontBuffer_.u[i * Ny_ + j + Nx_ * Ny_];
  velocity_outer[1]        = frontBuffer_.v[i * Ny_ + j + Nx_ * Ny_];
  velocity_outer[2]        = frontBuffer_.w[i * Ny_ + j + Nx_ * Ny_];
}
void Stencils::VectorBufferReadStencil::applyBackWall(int i, int j, int k) {
  RealType* velocity = vectorField_.getVector(i, j, k);
  velocity[0]        = backBuffer_.u[i * Ny_ + j];
  velocity[1]        = backBuffer_.v[i * Ny_ + j];
  velocity[2]        = backBuffer_.w[i * Ny_ + j];
}

RealType* Stencils::VectorBufferReadStencil::getLeftBuffer() { return leftBuffer_.u; }
RealType* Stencils::VectorBufferReadStencil::getRightBuffer() { return rightBuffer_.u; }
RealType* Stencils::VectorBufferReadStencil::getBottomBuffer() { return bottomBuffer_.u; }
RealType* Stencils::VectorBufferReadStencil::getTopBuffer() { return topBuffer_.u; }
RealType* Stencils::VectorBufferReadStencil::getFrontBuffer() { return frontBuffer_.u; }
RealType* Stencils::VectorBufferReadStencil::getBackBuffer() { return backBuffer_.u; }