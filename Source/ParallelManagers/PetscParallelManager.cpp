#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowField):
  Nx_(parameters.parallel.localSize[0] + 3),
  Ny_(parameters.parallel.localSize[1] + 3),
  Nz_(parameters.parallel.localSize[2] + 3),
  dim_(parameters.geometry.dim),
  leftNb_(parameters.parallel.leftNb),
  rightNb_(parameters.parallel.rightNb),
  bottomNb_(parameters.parallel.bottomNb),
  topNb_(parameters.parallel.topNb),
  frontNb_(parameters.parallel.frontNb),
  backNb_(parameters.parallel.backNb),
  pressureBufferFillStencil_(parameters, flowField.getPressure()),
  pressureBufferReadStencil_(parameters, flowField.getPressure()),
  velocityBufferFillStencil_(parameters, flowField.getVelocity()),
  velocityBufferReadStencil_(parameters, flowField.getVelocity()),
  parallelBoundaryIteratorPressureFill_(flowField, parameters, pressureBufferFillStencil_, 0, 0),
  parallelBoundaryIteratorPressureRead_(flowField, parameters, pressureBufferReadStencil_, 0, 0),
  parallelBoundaryIteratorVelocityFill_(flowField, parameters, velocityBufferFillStencil_, 0, 0),
  parallelBoundaryIteratorVelocityRead_(flowField, parameters, velocityBufferReadStencil_, 0, 0) {}


void ParallelManagers::PetscParallelManager::communicatePressure() {
  communicateBuffers(pressureBufferFillStencil_, pressureBufferReadStencil_, parallelBoundaryIteratorPressureFill_, parallelBoundaryIteratorPressureRead_, 0);
}
void ParallelManagers::PetscParallelManager::communicateVelocity() {
  communicateBuffers(velocityBufferFillStencil_, velocityBufferReadStencil_, parallelBoundaryIteratorVelocityFill_, parallelBoundaryIteratorVelocityRead_, 1);
}
void ParallelManagers::PetscParallelManager::communicateBuffers(
  Stencils::BufferStencil&             bufferFillStencil,
  Stencils::BufferStencil&             bufferReadStencil,
  ParallelBoundaryIterator<FlowField>& BoundaryIteratorFill,
  ParallelBoundaryIterator<FlowField>& BoundaryIteratorRead,
  bool                                 isVector
) {

  int buffMultiplier; // scalars have 1 buffer, vectors have dim buffers

  if (isVector == 1) {
    buffMultiplier = dim_;
  } else {
    buffMultiplier = 1;
  }

  if (dim_ == 2) {
    if (topNb_ >= 0) {
      BoundaryIteratorFill.iterateTop();
      MPI_Send(bufferFillStencil.getTopBuffer(), buffMultiplier * 2 * Nx_, MPI_DOUBLE, topNb_, 0, PETSC_COMM_WORLD);
    }
    if (bottomNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getBottomBuffer(), buffMultiplier * 2 * Nx_, MPI_DOUBLE, bottomNb_, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateBottom();
    }

    if (rightNb_ >= 0) {
      BoundaryIteratorFill.iterateRight();
      MPI_Send(bufferFillStencil.getRightBuffer(), buffMultiplier * 2 * Ny_, MPI_DOUBLE, rightNb_, 1, PETSC_COMM_WORLD);
    }
    if (leftNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getLeftBuffer(), buffMultiplier * 2 * Ny_, MPI_DOUBLE, leftNb_, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateLeft();
    }

    if (bottomNb_ >= 0) {
      BoundaryIteratorFill.iterateBottom();
      MPI_Send(bufferFillStencil.getBottomBuffer(), buffMultiplier * Nx_, MPI_DOUBLE, bottomNb_, 2, PETSC_COMM_WORLD);
    }
    if (topNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getTopBuffer(), buffMultiplier * Nx_, MPI_DOUBLE, topNb_, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateTop();
    }

    if (leftNb_ >= 0) {
      BoundaryIteratorFill.iterateLeft();
      MPI_Send(bufferFillStencil.getLeftBuffer(), buffMultiplier * Ny_, MPI_DOUBLE, leftNb_, 3, PETSC_COMM_WORLD);
    }
    if (rightNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getRightBuffer(), buffMultiplier * Ny_, MPI_DOUBLE, rightNb_, 3, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateRight();
    }
  } else if (dim_ == 3) {
    // Send top recive bottom
    if (topNb_ >= 0) {
      BoundaryIteratorFill.iterateTop();
      MPI_Send(bufferFillStencil.getTopBuffer(), buffMultiplier * 2 * Nx_ * Nz_, MPI_DOUBLE, topNb_, 0, PETSC_COMM_WORLD);
    }
    if (bottomNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getBottomBuffer(), buffMultiplier * 2 * Nz_ * Nx_, MPI_DOUBLE, bottomNb_, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateBottom();
    }
    // Send back recieve front
    if (backNb_ >= 0) {
      BoundaryIteratorFill.iterateBack();
      MPI_Send(bufferFillStencil.getBackBuffer(), buffMultiplier * 2 * Nx_ * Ny_, MPI_DOUBLE, backNb_, 1, PETSC_COMM_WORLD);
    }
    if (frontNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getFrontBuffer(), buffMultiplier * 2 * Nx_ * Ny_, MPI_DOUBLE, frontNb_, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateFront();
    }
    // Send right recieve left
    if (rightNb_ >= 0) {
      BoundaryIteratorFill.iterateRight();
      MPI_Send(bufferFillStencil.getRightBuffer(), buffMultiplier * 2 * Ny_ * Nz_, MPI_DOUBLE, rightNb_, 2, PETSC_COMM_WORLD);
    }
    if (leftNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getLeftBuffer(), buffMultiplier * 2 * Ny_ * Nz_, MPI_DOUBLE, leftNb_, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateLeft();
    }
    // Send bottom recieve top
    if (bottomNb_ >= 0) {
      BoundaryIteratorFill.iterateBottom();
      MPI_Send(bufferFillStencil.getBottomBuffer(), buffMultiplier * Nx_ * Nz_, MPI_DOUBLE, bottomNb_, 3, PETSC_COMM_WORLD);
    }
    if (topNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getTopBuffer(), buffMultiplier * Nx_ * Nz_, MPI_DOUBLE, topNb_, 3, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateTop();
    }
    // Send front recieve back
    if (frontNb_ >= 0) {
      BoundaryIteratorFill.iterateFront();
      MPI_Send(bufferFillStencil.getFrontBuffer(), buffMultiplier * Nx_ * Ny_, MPI_DOUBLE, frontNb_, 4, PETSC_COMM_WORLD);
    }
    if (backNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getBackBuffer(), buffMultiplier * Nx_ * Ny_, MPI_DOUBLE, backNb_, 4, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateBack();
    }
    // Send left recieve right
    if (leftNb_ >= 0) {
      BoundaryIteratorFill.iterateLeft();
      MPI_Send(bufferFillStencil.getLeftBuffer(), buffMultiplier * Ny_ * Nz_, MPI_DOUBLE, leftNb_, 5, PETSC_COMM_WORLD);
    }
    if (rightNb_ >= 0) {
      MPI_Recv(bufferReadStencil.getRightBuffer(), buffMultiplier * Ny_ * Nz_, MPI_DOUBLE, rightNb_, 5, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      BoundaryIteratorRead.iterateRight();
    }
  }
}
