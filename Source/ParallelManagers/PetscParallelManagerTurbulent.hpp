#pragma once

#include "ParallelManagers/PetscParallelManager.hpp"

// Class that builds on PetscParallelManager to add communication of the turbulent viscosity
namespace ParallelManagers {
  class PetscParallelManagerTurbulent: public PetscParallelManager {
  private:
    Stencils::ScalarBufferFillStencil   turbViscosityBufferFillStencil_;
    Stencils::ScalarBufferReadStencil   turbViscosityBufferReadStencil_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorTurbViscosityFill_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorTurbViscosityRead_;

    Stencils::ScalarBufferFillStencil   turbKineticEnergyBufferFillStencil_;
    Stencils::ScalarBufferReadStencil   turbKineticEnergyBufferReadStencil_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorTurbKineticEnergyFill_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorTurbKineticEnergyRead_;

    Stencils::ScalarBufferFillStencil   turbOmegaBufferFillStencil_;
    Stencils::ScalarBufferReadStencil   turbOmegaBufferReadStencil_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorTurbOmegaFill_;
    ParallelBoundaryIterator<FlowField> parallelBoundaryIteratorTurbOmegaRead_;

  public:
    PetscParallelManagerTurbulent(Parameters& parameters, TurbulentFlowField& flowField);
    ~PetscParallelManagerTurbulent() = default;

    void communicateTurbViscosity();
    void communicateTurbKineticEnergy();
    void communicateTurbOmega();
  };
} // namespace ParallelManagers