#include "StdAfx.hpp"

#include "PetscParallelManagerTurbulent.hpp"

ParallelManagers::PetscParallelManagerTurbulent::PetscParallelManagerTurbulent(Parameters& parameters, TurbulentFlowField& flowField):
  PetscParallelManager(parameters, flowField),
  turbViscosityBufferFillStencil_(parameters, flowField.getTurbViscosity()),
  turbViscosityBufferReadStencil_(parameters, flowField.getTurbViscosity()),
  parallelBoundaryIteratorTurbViscosityFill_(flowField, parameters, turbViscosityBufferFillStencil_, 0, 0),
  parallelBoundaryIteratorTurbViscosityRead_(flowField, parameters, turbViscosityBufferReadStencil_, 0, 0),
  turbKineticEnergyBufferFillStencil_(parameters, flowField.getTurbKineticEnergy()),
  turbKineticEnergyBufferReadStencil_(parameters, flowField.getTurbKineticEnergy()),
  parallelBoundaryIteratorTurbKineticEnergyFill_(flowField, parameters, turbKineticEnergyBufferFillStencil_, 0, 0),
  parallelBoundaryIteratorTurbKineticEnergyRead_(flowField, parameters, turbKineticEnergyBufferReadStencil_, 0, 0),
  turbOmegaBufferFillStencil_(parameters, flowField.getTurbOmega()),
  turbOmegaBufferReadStencil_(parameters, flowField.getTurbOmega()),
  parallelBoundaryIteratorTurbOmegaFill_(flowField, parameters, turbOmegaBufferFillStencil_, 0, 0),
  parallelBoundaryIteratorTurbOmegaRead_(flowField, parameters, turbOmegaBufferReadStencil_, 0, 0) {}

void ParallelManagers::PetscParallelManagerTurbulent::communicateTurbViscosity() {
  communicateBuffers(turbViscosityBufferFillStencil_, turbViscosityBufferReadStencil_, parallelBoundaryIteratorTurbViscosityFill_, parallelBoundaryIteratorTurbViscosityRead_, 0);
}

void ParallelManagers::PetscParallelManagerTurbulent::communicateTurbKineticEnergy() {
  communicateBuffers(
    turbKineticEnergyBufferFillStencil_, turbKineticEnergyBufferReadStencil_, parallelBoundaryIteratorTurbKineticEnergyFill_, parallelBoundaryIteratorTurbKineticEnergyRead_, 0
  );
}

void ParallelManagers::PetscParallelManagerTurbulent::communicateTurbOmega() {
  communicateBuffers(turbOmegaBufferFillStencil_, turbOmegaBufferReadStencil_, parallelBoundaryIteratorTurbOmegaFill_, parallelBoundaryIteratorTurbOmegaRead_, 0);
}
