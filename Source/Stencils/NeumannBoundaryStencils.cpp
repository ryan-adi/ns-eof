#include "StdAfx.hpp"

#include "NeumannBoundaryStencils.hpp"

Stencils::NeumannVelocityBoundaryStencil::NeumannVelocityBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

void Stencils::NeumannVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i - 1, j)[0] = flowField.getVelocity().getVector(i, j)[0];
  flowField.getVelocity().getVector(i, j)[1]     = flowField.getVelocity().getVector(i + 1, j)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i - 1, j)[0];
  flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0]     = flowField.getVelocity().getVector(i, j + 1)[0];
  flowField.getVelocity().getVector(i, j - 1)[1] = flowField.getVelocity().getVector(i, j)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i, j - 1)[0];
  flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i, j - 1)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i - 1, j, k)[0] = flowField.getVelocity().getVector(i, j, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i + 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i - 1, j, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i - 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j + 1, k)[0];
  flowField.getVelocity().getVector(i, j - 1, k)[1] = flowField.getVelocity().getVector(i, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j - 1, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j - 1, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j, k + 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i, j, k + 1)[1];
  flowField.getVelocity().getVector(i, j, k - 1)[2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j, k - 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j, k - 1)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j, k - 1)[2];
}

Stencils::NeumannFGHBoundaryStencil::NeumannFGHBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

// These are left empty. The right values should be computed by the FGH body stencil.
void Stencils::NeumannFGHBoundaryStencil::applyLeftWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::NeumannFGHBoundaryStencil::applyRightWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::NeumannFGHBoundaryStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::NeumannFGHBoundaryStencil::applyTopWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}

void Stencils::NeumannFGHBoundaryStencil::applyLeftWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::NeumannFGHBoundaryStencil::applyRightWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::NeumannFGHBoundaryStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::NeumannFGHBoundaryStencil::applyTopWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::NeumannFGHBoundaryStencil::applyFrontWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::NeumannFGHBoundaryStencil::applyBackWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}

// Turbulent Viscosity WS2 //////////////////////////////////////////////////////////////////
Stencils::NeumannTVBoundaryStencil::NeumannTVBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters) {}

// These are left empty. The right values should be computed by the FGH body stencil.
void Stencils::NeumannTVBoundaryStencil::applyLeftWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbViscosity().getScalar(i - 1, j) = flowField.getTurbViscosity().getScalar(i, j);
}
void Stencils::NeumannTVBoundaryStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbViscosity().getScalar(i, j) = flowField.getTurbViscosity().getScalar(i - 1, j);
}
void Stencils::NeumannTVBoundaryStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbViscosity().getScalar(i, j - 1) = flowField.getTurbViscosity().getScalar(i, j);
}
void Stencils::NeumannTVBoundaryStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbViscosity().getScalar(i, j) = flowField.getTurbViscosity().getScalar(i, j - 1);
}

void Stencils::NeumannTVBoundaryStencil::applyLeftWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbViscosity().getScalar(i - 1, j, k) = flowField.getTurbViscosity().getScalar(i, j, k);
}
void Stencils::NeumannTVBoundaryStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbViscosity().getScalar(i, j, k) = flowField.getTurbViscosity().getScalar(i - 1, j, k);
}
void Stencils::NeumannTVBoundaryStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbViscosity().getScalar(i, j - 1, k) = flowField.getTurbViscosity().getScalar(i, j, k);
}
void Stencils::NeumannTVBoundaryStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbViscosity().getScalar(i, j, k) = flowField.getTurbViscosity().getScalar(i, j - 1, k);
}
void Stencils::NeumannTVBoundaryStencil::applyFrontWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbViscosity().getScalar(i, j, k - 1) = flowField.getTurbViscosity().getScalar(i, j, k);
}
void Stencils::NeumannTVBoundaryStencil::applyBackWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbViscosity().getScalar(i, j, k) = flowField.getTurbViscosity().getScalar(i, j, k - 1);
}

// Turbulent kin energy WS3 //////////////////////////////////////////////////////////////////
Stencils::NeumannKBoundaryStencil::NeumannKBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters) {}

// These are left empty. The right values should be computed by the FGH body stencil.
void Stencils::NeumannKBoundaryStencil::applyLeftWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbKineticEnergy().getScalar(i - 1, j) = flowField.getTurbKineticEnergy().getScalar(i, j);
}
void Stencils::NeumannKBoundaryStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j) = flowField.getTurbKineticEnergy().getScalar(i - 1, j);
}
void Stencils::NeumannKBoundaryStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j - 1) = flowField.getTurbKineticEnergy().getScalar(i, j);
}
void Stencils::NeumannKBoundaryStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j) = flowField.getTurbKineticEnergy().getScalar(i, j - 1);
}

void Stencils::NeumannKBoundaryStencil::applyLeftWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbKineticEnergy().getScalar(i - 1, j, k) = flowField.getTurbKineticEnergy().getScalar(i, j, k);
}
void Stencils::NeumannKBoundaryStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbKineticEnergy().getScalar(i, j, k) = flowField.getTurbKineticEnergy().getScalar(i - 1, j, k);
}
void Stencils::NeumannKBoundaryStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbKineticEnergy().getScalar(i, j - 1, k) = flowField.getTurbKineticEnergy().getScalar(i, j, k);
}
void Stencils::NeumannKBoundaryStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbKineticEnergy().getScalar(i, j, k) = flowField.getTurbKineticEnergy().getScalar(i, j - 1, k);
}
void Stencils::NeumannKBoundaryStencil::applyFrontWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbKineticEnergy().getScalar(i, j, k - 1) = flowField.getTurbKineticEnergy().getScalar(i, j, k);
}
void Stencils::NeumannKBoundaryStencil::applyBackWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbKineticEnergy().getScalar(i, j, k) = flowField.getTurbKineticEnergy().getScalar(i, j, k - 1);
}

// Turbulent omega WS3 //////////////////////////////////////////////////////////////////
Stencils::NeumannWBoundaryStencil::NeumannWBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters) {}

// These are left empty. The right values should be computed by the FGH body stencil.
void Stencils::NeumannWBoundaryStencil::applyLeftWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbOmega().getScalar(i - 1, j) = flowField.getTurbOmega().getScalar(i, j);
}
void Stencils::NeumannWBoundaryStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbOmega().getScalar(i, j) = flowField.getTurbOmega().getScalar(i - 1, j);
}
void Stencils::NeumannWBoundaryStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbOmega().getScalar(i, j - 1) = flowField.getTurbOmega().getScalar(i, j);
}
void Stencils::NeumannWBoundaryStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {
  flowField.getTurbOmega().getScalar(i, j) = flowField.getTurbOmega().getScalar(i, j - 1);
}

void Stencils::NeumannWBoundaryStencil::applyLeftWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbOmega().getScalar(i - 1, j, k) = flowField.getTurbOmega().getScalar(i, j, k);
}
void Stencils::NeumannWBoundaryStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = flowField.getTurbOmega().getScalar(i - 1, j, k);
}
void Stencils::NeumannWBoundaryStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbOmega().getScalar(i, j - 1, k) = flowField.getTurbOmega().getScalar(i, j, k);
}
void Stencils::NeumannWBoundaryStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = flowField.getTurbOmega().getScalar(i, j - 1, k);
}
void Stencils::NeumannWBoundaryStencil::applyFrontWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbOmega().getScalar(i, j, k - 1) = flowField.getTurbOmega().getScalar(i, j, k);
}
void Stencils::NeumannWBoundaryStencil::applyBackWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = flowField.getTurbOmega().getScalar(i, j, k - 1);
}