#pragma once

#include "Definitions.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"
#include "Stencils/ScalarBufferFillStencil.hpp"
#include "Stencils/ScalarBufferReadStencil.hpp"
#include "Stencils/VectorBufferFillStencil.hpp"
#include "Stencils/VectorBufferReadStencil.hpp"

namespace ParallelManagers {
  class PetscParallelManager {
  private:
    int Nx_, Ny_, Nz_;
    int dim_      = -1;
    int leftNb_   = MPI_PROC_NULL;
    int rightNb_  = MPI_PROC_NULL;
    int bottomNb_ = MPI_PROC_NULL;
    int topNb_    = MPI_PROC_NULL;
    int frontNb_  = MPI_PROC_NULL;
    int backNb_   = MPI_PROC_NULL;

    Stencils::ScalarBufferFillStencil   pressureBufferFillStencil_;
    Stencils::ScalarBufferReadStencil   pressureBufferReadStencil_;
    Stencils::VectorBufferFillStencil   velocityBufferFillStencil_;
    Stencils::VectorBufferReadStencil   velocityBufferReadStencil_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorPressureFill_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorPressureRead_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorVelocityFill_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorVelocityRead_;

  public:
    PetscParallelManager(Parameters& parameters, FlowField& flowField);
    ~PetscParallelManager() = default;

    void communicateBuffers(
      Stencils::BufferStencil&             bufferFillStencil,
      Stencils::BufferStencil&             bufferReadStencil,
      ParallelBoundaryIterator<FlowField>& BoundaryIteratorFill,
      ParallelBoundaryIterator<FlowField>& BoundaryIteratorRead,
      bool                                 isVector
    );
    void communicatePressure();
    void communicateVelocity();
  };
} // namespace ParallelManagers