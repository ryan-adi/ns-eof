#include "StdAfx.hpp"

#include "VectorBufferFillStencil.hpp"

Stencils::VectorBufferFillStencil::VectorBufferFillStencil(const Parameters& parameters, VectorField& vectorField):
  BufferStencil(parameters),
  vectorField_(vectorField),
  leftBuffer_((parameters.geometry.dim == 3 ? Ny_ * Nz_ : Ny_), parameters.geometry.dim),
  rightBuffer_(2 * (parameters.geometry.dim == 3 ? Ny_ * Nz_ : Ny_), parameters.geometry.dim),
  bottomBuffer_((parameters.geometry.dim == 3 ? Nx_ * Nz_ : Nx_), parameters.geometry.dim),
  topBuffer_(2 * (parameters.geometry.dim == 3 ? Nx_ * Nz_ : Nx_), parameters.geometry.dim),
  frontBuffer_((parameters.geometry.dim == 3 ? Nx_ * Ny_ : 0), parameters.geometry.dim),
  backBuffer_(2 * (parameters.geometry.dim == 3 ? Nx_ * Ny_ : 0), parameters.geometry.dim) {}

Stencils::VectorBufferFillStencil::~VectorBufferFillStencil() {}
// example double buffer, for right wall
// 1 2
// 1 2
// 1 2
// 1 2
// velocity has the 2 values, velocity_left has the 1 values
// in the buffer they are stored as 2 2 2 2 1 1 1 1
// so when we read it from the buffer of the right neighbour we have to copy
// the 2 as the inner values and the 1 as the outer values. I.e. we copy the first
// Ny_ elements from the buffer to the inner values and the last Ny_ elements to the
// outer values

// example double buffer, for top wall
// 2 2 2 2
// 1 1 1 1
// velocity has the 2 values, velocity_bottom has the 1 values
// in the buffer they are stored as 2 2 2 2 1 1 1 1
// so when we read it from the buffer of the top neighbour we have to copy
// the 2 as the inner values and the 1 as the outer values. I.e. we copy the first
// Nx_ elements from the buffer to the inner values and the last Nx_ elements to the
// outer values

void Stencils::VectorBufferFillStencil::applyLeftWall(int i, int j) {
  RealType* velocity = vectorField_.getVector(i + 2, j);
  leftBuffer_.u[j]   = velocity[0];
  leftBuffer_.v[j]   = velocity[1];
}
void Stencils::VectorBufferFillStencil::applyRightWall(int i, int j) {
  RealType* velocity      = vectorField_.getVector(i - 1, j);
  rightBuffer_.u[j]       = velocity[0];
  rightBuffer_.v[j]       = velocity[1];
  RealType* velocity_left = vectorField_.getVector(i - 2, j);
  rightBuffer_.u[j + Ny_] = velocity_left[0];
  rightBuffer_.v[j + Ny_] = velocity_left[1];
}
void Stencils::VectorBufferFillStencil::applyBottomWall(int i, int j) {
  RealType* velocity = vectorField_.getVector(i, j + 2);
  bottomBuffer_.u[i] = velocity[0];
  bottomBuffer_.v[i] = velocity[1];
}
void Stencils::VectorBufferFillStencil::applyTopWall(int i, int j) {
  RealType* velocity        = vectorField_.getVector(i, j - 1);
  topBuffer_.u[i]           = velocity[0];
  topBuffer_.v[i]           = velocity[1];
  RealType* velocity_bottom = vectorField_.getVector(i, j - 2);
  topBuffer_.u[i + Nx_]     = velocity_bottom[0];
  topBuffer_.v[i + Nx_]     = velocity_bottom[1];
}

void Stencils::VectorBufferFillStencil::applyLeftWall(int i, int j, int k) {
  RealType* velocity         = vectorField_.getVector(i + 2, j, k);
  leftBuffer_.u[j * Nz_ + k] = velocity[0];
  leftBuffer_.v[j * Nz_ + k] = velocity[1];
  leftBuffer_.w[j * Nz_ + k] = velocity[2];
}
void Stencils::VectorBufferFillStencil::applyRightWall(int i, int j, int k) {
  RealType* velocity                      = vectorField_.getVector(i - 1, j, k);
  rightBuffer_.u[j * Nz_ + k]             = velocity[0];
  rightBuffer_.v[j * Nz_ + k]             = velocity[1];
  rightBuffer_.w[j * Nz_ + k]             = velocity[2];
  RealType* velocity_left                 = vectorField_.getVector(i - 2, j, k);
  rightBuffer_.u[j * Nz_ + k + Ny_ * Nz_] = velocity_left[0];
  rightBuffer_.v[j * Nz_ + k + Ny_ * Nz_] = velocity_left[1];
  rightBuffer_.w[j * Nz_ + k + Ny_ * Nz_] = velocity_left[2];
}
void Stencils::VectorBufferFillStencil::applyBottomWall(int i, int j, int k) {
  RealType* velocity           = vectorField_.getVector(i, j + 2, k);
  bottomBuffer_.u[i * Nz_ + k] = velocity[0];
  bottomBuffer_.v[i * Nz_ + k] = velocity[1];
  bottomBuffer_.w[i * Nz_ + k] = velocity[2];
}
void Stencils::VectorBufferFillStencil::applyTopWall(int i, int j, int k) {
  RealType* velocity                    = vectorField_.getVector(i, j - 1, k);
  topBuffer_.u[i * Nz_ + k]             = velocity[0];
  topBuffer_.v[i * Nz_ + k]             = velocity[1];
  topBuffer_.w[i * Nz_ + k]             = velocity[2];
  RealType* velocity_bottom             = vectorField_.getVector(i, j - 2, k);
  topBuffer_.u[i * Nz_ + k + Nx_ * Nz_] = velocity_bottom[0];
  topBuffer_.v[i * Nz_ + k + Nx_ * Nz_] = velocity_bottom[1];
  topBuffer_.w[i * Nz_ + k + Nx_ * Nz_] = velocity_bottom[2];
}
void Stencils::VectorBufferFillStencil::applyFrontWall(int i, int j, int k) {
  RealType* velocity          = vectorField_.getVector(i, j, k + 2);
  frontBuffer_.u[i * Ny_ + j] = velocity[0];
  frontBuffer_.v[i * Ny_ + j] = velocity[1];
  frontBuffer_.w[i * Ny_ + j] = velocity[2];
}
void Stencils::VectorBufferFillStencil::applyBackWall(int i, int j, int k) {
  RealType* velocity                     = vectorField_.getVector(i, j, k - 1);
  backBuffer_.u[i * Ny_ + j]             = velocity[0];
  backBuffer_.v[i * Ny_ + j]             = velocity[1];
  backBuffer_.w[i * Ny_ + j]             = velocity[2];
  RealType* velocity_front               = vectorField_.getVector(i, j, k - 2);
  backBuffer_.u[i * Ny_ + j + Nx_ * Ny_] = velocity_front[0];
  backBuffer_.v[i * Ny_ + j + Nx_ * Ny_] = velocity_front[1];
  backBuffer_.w[i * Ny_ + j + Nx_ * Ny_] = velocity_front[2];
}

RealType* Stencils::VectorBufferFillStencil::getLeftBuffer() { return leftBuffer_.u; }
RealType* Stencils::VectorBufferFillStencil::getRightBuffer() { return rightBuffer_.u; }
RealType* Stencils::VectorBufferFillStencil::getBottomBuffer() { return bottomBuffer_.u; }
RealType* Stencils::VectorBufferFillStencil::getTopBuffer() { return topBuffer_.u; }
RealType* Stencils::VectorBufferFillStencil::getFrontBuffer() { return frontBuffer_.u; }
RealType* Stencils::VectorBufferFillStencil::getBackBuffer() { return backBuffer_.u; }
