#include "StdAfx.hpp"

#include "MovingWallStencils.hpp"

Stencils::MovingWallVelocityStencil::MovingWallVelocityStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

void Stencils::MovingWallVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = parameters_.walls.vectorLeft[0];
  flowField.getVelocity().getVector(i, j)[1] = 2 * parameters_.walls.vectorLeft[1] - flowField.getVelocity().getVector(i + 1, j)[1];
}

void Stencils::MovingWallVelocityStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i - 1, j)[0] = parameters_.walls.vectorRight[0];
  flowField.getVelocity().getVector(i, j)[1]     = 2 * parameters_.walls.vectorRight[1] - flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::MovingWallVelocityStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = 2 * parameters_.walls.vectorBottom[0] - flowField.getVelocity().getVector(i, j + 1)[0];
  flowField.getVelocity().getVector(i, j)[1] = parameters_.walls.vectorBottom[1];
}

void Stencils::MovingWallVelocityStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0]     = 2 * parameters_.walls.vectorTop[0] - flowField.getVelocity().getVector(i, j - 1)[0];
  flowField.getVelocity().getVector(i, j - 1)[1] = parameters_.walls.vectorTop[1];
}

void Stencils::MovingWallVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = parameters_.walls.vectorLeft[0];
  flowField.getVelocity().getVector(i, j, k)[1] = 2 * parameters_.walls.vectorLeft[1] - flowField.getVelocity().getVector(i + 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = 2 * parameters_.walls.vectorLeft[2] - flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i - 1, j, k)[0] = parameters_.walls.vectorRight[0];
  flowField.getVelocity().getVector(i, j, k)[1]     = 2 * parameters_.walls.vectorRight[1] - flowField.getVelocity().getVector(i - 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2]     = 2 * parameters_.walls.vectorRight[2] - flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = 2 * parameters_.walls.vectorBottom[0] - flowField.getVelocity().getVector(i, j + 1, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = parameters_.walls.vectorBottom[1];
  flowField.getVelocity().getVector(i, j, k)[2] = 2 * parameters_.walls.vectorBottom[2] - flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = 2 * parameters_.walls.vectorTop[0] - flowField.getVelocity().getVector(i, j - 1, k)[0];
  flowField.getVelocity().getVector(i, j - 1, k)[1] = parameters_.walls.vectorTop[1];
  flowField.getVelocity().getVector(i, j, k)[2]     = 2 * parameters_.walls.vectorTop[2] - flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = 2 * parameters_.walls.vectorFront[0] - flowField.getVelocity().getVector(i, j, k + 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = 2 * parameters_.walls.vectorFront[1] - flowField.getVelocity().getVector(i, j, k + 1)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = parameters_.walls.vectorFront[2];
}

void Stencils::MovingWallVelocityStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = 2 * parameters_.walls.vectorBack[0] - flowField.getVelocity().getVector(i, j, k - 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1]     = 2 * parameters_.walls.vectorBack[1] - flowField.getVelocity().getVector(i, j, k - 1)[1];
  flowField.getVelocity().getVector(i, j, k - 1)[2] = parameters_.walls.vectorBack[2];
}

Stencils::MovingWallFGHStencil::MovingWallFGHStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

void Stencils::MovingWallFGHStencil::applyLeftWall(FlowField& flowField, int i, int j) { flowField.getFGH().getVector(i, j)[0] = parameters_.walls.vectorLeft[0]; }

void Stencils::MovingWallFGHStencil::applyRightWall(FlowField& flowField, int i, int j) { flowField.getFGH().getVector(i - 1, j)[0] = parameters_.walls.vectorRight[0]; }

void Stencils::MovingWallFGHStencil::applyBottomWall(FlowField& flowField, int i, int j) { flowField.getFGH().getVector(i, j)[1] = parameters_.walls.vectorBottom[1]; }

void Stencils::MovingWallFGHStencil::applyTopWall(FlowField& flowField, int i, int j) { flowField.getFGH().getVector(i, j - 1)[1] = parameters_.walls.vectorTop[1]; }

void Stencils::MovingWallFGHStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) { flowField.getFGH().getVector(i, j, k)[0] = parameters_.walls.vectorLeft[0]; }

void Stencils::MovingWallFGHStencil::applyRightWall(FlowField& flowField, int i, int j, int k) { flowField.getFGH().getVector(i - 1, j, k)[0] = parameters_.walls.vectorRight[0]; }

void Stencils::MovingWallFGHStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) { flowField.getFGH().getVector(i, j, k)[1] = parameters_.walls.vectorBottom[1]; }

void Stencils::MovingWallFGHStencil::applyTopWall(FlowField& flowField, int i, int j, int k) { flowField.getFGH().getVector(i, j - 1, k)[1] = parameters_.walls.vectorTop[1]; }

void Stencils::MovingWallFGHStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) { flowField.getFGH().getVector(i, j, k)[2] = parameters_.walls.vectorFront[2]; }

void Stencils::MovingWallFGHStencil::applyBackWall(FlowField& flowField, int i, int j, int k) { flowField.getFGH().getVector(i, j, k - 1)[2] = parameters_.walls.vectorBack[2]; }

// Turbulent Viscosity WS2 //////////////////////////////////////////////////////////////////
Stencils::MovingWallTVStencil::MovingWallTVStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters) {}

void Stencils::MovingWallTVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbViscosity().getScalar(i, j) = -flowField.getTurbViscosity().getScalar(i + 1, j);
}
void Stencils::MovingWallTVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbViscosity().getScalar(i, j) = -flowField.getTurbViscosity().getScalar(i - 1, j);
}
void Stencils::MovingWallTVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbViscosity().getScalar(i, j) = -flowField.getTurbViscosity().getScalar(i, j + 1);
}
void Stencils::MovingWallTVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbViscosity().getScalar(i, j) = -flowField.getTurbViscosity().getScalar(i, j - 1);
}
void Stencils::MovingWallTVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbViscosity().getScalar(i, j, k) = 0.0; }
void Stencils::MovingWallTVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbViscosity().getScalar(i, j, k) = 0.0; }
void Stencils::MovingWallTVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbViscosity().getScalar(i, j, k) = 0.0; }
void Stencils::MovingWallTVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbViscosity().getScalar(i, j, k) = 0.0; }
void Stencils::MovingWallTVStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbViscosity().getScalar(i, j, k) = 0.0; }
void Stencils::MovingWallTVStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbViscosity().getScalar(i, j, k) = 0.0; }

//  Turb Kinetic Energy K WS3 //////////////////////////////////////////////////////////////////
Stencils::MovingWallKStencil::MovingWallKStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters) {}

void Stencils::MovingWallKStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j) = -flowField.getTurbKineticEnergy().getScalar(i + 1, j);
}
void Stencils::MovingWallKStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j) = -flowField.getTurbKineticEnergy().getScalar(i - 1, j);
}
void Stencils::MovingWallKStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j) = -flowField.getTurbKineticEnergy().getScalar(i, j + 1);
}
void Stencils::MovingWallKStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j) = -flowField.getTurbKineticEnergy().getScalar(i, j - 1);
}
void Stencils::MovingWallKStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbKineticEnergy().getScalar(i, j, k) = 0; }
void Stencils::MovingWallKStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbKineticEnergy().getScalar(i, j, k) = 0; }
void Stencils::MovingWallKStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbKineticEnergy().getScalar(i, j, k) = 0; }
void Stencils::MovingWallKStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbKineticEnergy().getScalar(i, j, k) = 0; }
void Stencils::MovingWallKStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbKineticEnergy().getScalar(i, j, k) = 0; }
void Stencils::MovingWallKStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) { flowField.getTurbKineticEnergy().getScalar(i, j, k) = 0; }

//  Turb Omega W WS3 //////////////////////////////////////////////////////////////////
RealType computeW2D(RealType Re, RealType Dx) {
  RealType beta1 = 3.0 / 40.0;
  return 6.0 / (Re * beta1 * Dx * Dx);
}

Stencils::MovingWallWStencil::MovingWallWStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters) {}

void Stencils::MovingWallWStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbOmega().getScalar(i, j) = flowField.getTurbOmega().getScalar(i + 1, j);
}

void Stencils::MovingWallWStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbOmega().getScalar(i, j) = flowField.getTurbOmega().getScalar(i - 1, j);
}

void Stencils::MovingWallWStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  RealType Dy      = parameters_.meshsize->getDy(i, j);
  RealType Dy_1    = parameters_.meshsize->getDy(i, j + 1);
  RealType Re      = parameters_.flow.Re;
  RealType omega_w = computeW2D(Re, Dy_1 / 2);
  RealType omega_1 = flowField.getTurbOmega().getScalar(i, j + 1);
  flowField.getTurbOmega().getScalar(i, j) = ((Dy + Dy_1) / Dy_1) * (omega_w - omega_1 * (Dy / (Dy + Dy_1)));
}

void Stencils::MovingWallWStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  RealType Dy      = parameters_.meshsize->getDy(i, j);
  RealType Dy_1    = parameters_.meshsize->getDy(i, j - 1);
  RealType Re      = parameters_.flow.Re;
  RealType omega_w = computeW2D(Re, Dy_1 / 2);
  RealType omega_1 = flowField.getTurbOmega().getScalar(i, j - 1);
  flowField.getTurbOmega().getScalar(i, j) = ((Dy + Dy_1) / Dy_1) * (omega_w - omega_1 * (Dy / (Dy + Dy_1)));
}

void Stencils::MovingWallWStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = flowField.getTurbOmega().getScalar(i + 1, j, k);
}

void Stencils::MovingWallWStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = flowField.getTurbOmega().getScalar(i - 1, j, k);
}

void Stencils::MovingWallWStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType Dy                              = parameters_.meshsize->getDy(i, j, k);
  RealType Dy_1                            = parameters_.meshsize->getDy(i, j + 1, k);
  RealType Re                              = parameters_.flow.Re;
  RealType beta1                           = 3.0 / 40.0;
  flowField.getTurbOmega().getScalar(i, j) = 6.0 / (Re * beta1 * (Dy + Dy_1) * (Dy + Dy_1) * 0.25);
}

void Stencils::MovingWallWStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType Dy                              = parameters_.meshsize->getDy(i, j, k);
  RealType Dy_1                            = parameters_.meshsize->getDy(i, j + 1, k);
  RealType Re                              = parameters_.flow.Re;
  RealType beta1                           = 3.0 / 40.0;
  flowField.getTurbOmega().getScalar(i, j) = 6.0 / (Re * beta1 * (Dy + Dy_1) * (Dy + Dy_1) * 0.25);
}

void Stencils::MovingWallWStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = flowField.getTurbOmega().getScalar(i, j, k + 1);
}

void Stencils::MovingWallWStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = flowField.getTurbOmega().getScalar(i, j, k - 1);
}