#include "StdAfx.hpp"

#include "BFInputStencils.hpp"

RealType computeVelocity3D(int i, int j, int k, RealType stepSize, const Parameters& parameters) {
  const RealType posY = parameters.meshsize->getPosY(i, j, k);
  const RealType posZ = parameters.meshsize->getPosZ(i, j, k);
  const RealType dy   = parameters.meshsize->getDy(i, j, k);
  const RealType dz   = parameters.meshsize->getDz(i, j, k);

  if (posY + 0.5 * dy >= stepSize) {
    // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells.
    // const RealType inletYSize = parameters.geometry.lengthY - stepSize;
    // const RealType inletZSize = parameters.geometry.lengthZ;

    // const RealType y = posY + 0.5 * dy - stepSize;
    // const RealType z = posZ + 0.5 * dz;

    // return 36.0 * parameters.walls.vectorLeft[0] / (inletZSize * inletZSize * inletYSize * inletYSize) * y * (y - inletYSize) * z * (z - inletZSize);
    // block profile
    return parameters.walls.vectorLeft[0];
  } else {
    return 0.0;
  }
}

RealType computeVelocity2D(int i, int j, RealType stepSize, const Parameters& parameters) {
  const RealType posY = parameters.meshsize->getPosY(i, j);
  const RealType dy   = parameters.meshsize->getDy(i, j);


  if (posY + 0.5 * dy >= stepSize) {
    // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells.
    // const RealType inletYSize = parameters.geometry.lengthY - stepSize;

    // const RealType y = posY + 0.5 * dy - stepSize;

    // For turbulence, please use: return parameters.walls.vectorLeft[0]; WS2
    // parabolic profile
    // return 6.0 * parameters.walls.vectorLeft[0] / (inletYSize * inletYSize) * y * (inletYSize - y);

    // block profile
    return parameters.walls.vectorLeft[0];

  } else {
    return 0.0;
  }
}

RealType computeTKE2D([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] RealType stepSize, [[maybe_unused]] const Parameters& parameters) {
  const RealType Re_dh         = parameters.flow.Re * (1.0 - parameters.bfStep.yRatio) / parameters.geometry.lengthY;
  const RealType turbIntensity = 0.16 * std::pow(Re_dh, -0.125);
  const RealType turbKinEnergy = 3.0 / 2.0 * turbIntensity * turbIntensity; // 3/2 * (U*I)^2 with U = 1.0
  return turbKinEnergy;
}
RealType computeTKE3D([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k, [[maybe_unused]] RealType stepSize, [[maybe_unused]] const Parameters& parameters) {
  const RealType Re_dh         = parameters.flow.Re * (1.0 - parameters.bfStep.yRatio) / parameters.geometry.lengthY;
  const RealType turbIntensity = 0.16 * std::pow(Re_dh, -0.125);
  const RealType turbKinEnergy = 3.0 / 2.0 * turbIntensity * turbIntensity; // 3/2 * (U*I)^2 with U = 1.0
  return turbKinEnergy;
}

RealType computeW2D([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] RealType stepSize, const Parameters& parameters) {
  const RealType Re_dh         = parameters.flow.Re * (1.0 - parameters.bfStep.yRatio) / parameters.geometry.lengthY;
  const RealType turbIntensity = 0.16 * std::pow(Re_dh, -0.125);
  const RealType turbLenght    = 0.07 * parameters.geometry.lengthY;                           // 7% of the hydraulic diameter, for channel flow it is the height of the channel
  const RealType turbKinEnergy = 3.0 / 2.0 * turbIntensity * turbIntensity;                    // 3/2 * (U*I)^2 with U = 1.0
  const RealType turbOmega     = std::pow(0.09, -0.25) * std::sqrt(turbKinEnergy) / turbLenght; // C_mu^1/4 * sqrt(k) / l
  return turbOmega;
}
RealType computeW3D([[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k, [[maybe_unused]] RealType stepSize, const Parameters& parameters) {
  const RealType Re_dh         = parameters.flow.Re * (1.0 - parameters.bfStep.yRatio) / parameters.geometry.lengthY;
  const RealType turbIntensity = 0.16 * std::pow(Re_dh, -0.125);
  const RealType turbLenght    = 0.07 * parameters.geometry.lengthY;                           // 7% of the hydraulic diameter, for channel flow it is the height of the channel
  const RealType turbKinEnergy = 3.0 / 2.0 * turbIntensity * turbIntensity;                    // 3/2 * (U*I)^2 with U = 1.0
  const RealType turbOmega     = std::pow(0.09, -0.25) * std::sqrt(turbKinEnergy) / turbLenght; // C_mu^1/4 * sqrt(k) / l
  return turbOmega;
}


Stencils::BFInputVelocityStencil::BFInputVelocityStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters)
  // Here, the obstacle size is set to zero if it was set as negative at the configuration.
  ,
  stepSize_(parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0) {

  if (parameters_.geometry.dim == 2) {
    RealType posY   = parameters_.meshsize->getPosY(0, 0);
    RealType dy     = parameters_.meshsize->getDy(0, 0);
    RealType nextDy = parameters_.meshsize->getDy(0, 1);

    for (int j = 0; j < parameters_.geometry.sizeY - 1; ++j) {
      posY   = parameters_.meshsize->getPosY(0, j);
      dy     = parameters_.meshsize->getDy(0, j);
      nextDy = parameters_.meshsize->getDy(0, j + 1);

      // Check if stepSize is in this cell
      if (posY + 0.5 * dy < stepSize_ && stepSize_ <= posY + dy + 0.5 * nextDy) {
        stepSize_ = posY + dy;
        break;
      }
    }
  } else if (parameters_.geometry.dim == 3) {
    RealType posY   = parameters_.meshsize->getPosY(0, 0, 0);
    RealType dy     = parameters_.meshsize->getDy(0, 0, 0);
    RealType nextDy = parameters_.meshsize->getDy(0, 1, 0);

    for (int j = 0; j < parameters_.geometry.sizeY - 1; ++j) {
      posY   = parameters_.meshsize->getPosY(0, j, 0);
      dy     = parameters_.meshsize->getDy(0, j, 0);
      nextDy = parameters_.meshsize->getDy(0, j + 1, 0);

      if (posY + 0.5 * dy < stepSize_ && stepSize_ <= posY + dy + 0.5 * nextDy) {
        stepSize_ = posY + dy;
        break;
      }
    }
  }
}

void Stencils::BFInputVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = computeVelocity2D(i, j, stepSize_, parameters_);
  flowField.getVelocity().getVector(i, j)[1] = -flowField.getVelocity().getVector(i + 1, j)[1];
}

// Most of the functions are empty, and they shouldn't be called, assuming that the input is always located at the left.
void Stencils::BFInputVelocityStencil::applyRightWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputVelocityStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputVelocityStencil::applyTopWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}

void Stencils::BFInputVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = computeVelocity3D(i, j, k, stepSize_, parameters_);
  flowField.getVelocity().getVector(i, j, k)[1] = -flowField.getVelocity().getVector(i + 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = -flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::BFInputVelocityStencil::applyRightWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputVelocityStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputVelocityStencil::applyTopWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputVelocityStencil::applyFrontWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputVelocityStencil::applyBackWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}

Stencils::BFInputFGHStencil::BFInputFGHStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters),
  stepSize_(parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0) {}

void Stencils::BFInputFGHStencil::applyLeftWall(FlowField& flowField, int i, int j) { flowField.getFGH().getVector(i, j)[0] = computeVelocity2D(i, j, stepSize_, parameters_); }

void Stencils::BFInputFGHStencil::applyRightWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputFGHStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputFGHStencil::applyTopWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}

void Stencils::BFInputFGHStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i, j, k)[0] = computeVelocity3D(i, j, k, stepSize_, parameters_);
}

void Stencils::BFInputFGHStencil::applyRightWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputFGHStencil::applyBottomWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputFGHStencil::applyTopWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputFGHStencil::applyFrontWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
void Stencils::BFInputFGHStencil::applyBackWall([[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}

// WS2 Turbulent Viscosity input, does nothing
Stencils::BFInputTVStencil::BFInputTVStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  stepSize_(parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0) {}

void Stencils::BFInputTVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbViscosity().getScalar(i, j) = flowField.getTurbKineticEnergy().getScalar(i, j) / flowField.getTurbOmega().getScalar(i, j);
}

void Stencils::BFInputTVStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputTVStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputTVStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}

void Stencils::BFInputTVStencil::applyLeftWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}

void Stencils::BFInputTVStencil::applyRightWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputTVStencil::applyBottomWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputTVStencil::applyTopWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputTVStencil::applyFrontWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputTVStencil::applyBackWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// WS3 Turbulent Kinetic Energy input
Stencils::BFInputKStencil::BFInputKStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  stepSize_(parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0) {}

void Stencils::BFInputKStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  flowField.getTurbKineticEnergy().getScalar(i, j) = computeTKE2D(i, j, stepSize_, parameters_);
}

void Stencils::BFInputKStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputKStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputKStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}

void Stencils::BFInputKStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getTurbKineticEnergy().getScalar(i, j, k) = computeTKE3D(i, j, k, stepSize_, parameters_);
}

void Stencils::BFInputKStencil::applyRightWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputKStencil::applyBottomWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputKStencil::applyTopWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputKStencil::applyFrontWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputKStencil::applyBackWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// WS3 Turbulent omega input
Stencils::BFInputWStencil::BFInputWStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  stepSize_(parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0) {}

void Stencils::BFInputWStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) { flowField.getTurbOmega().getScalar(i, j) = computeW2D(i, j, stepSize_, parameters_); }

void Stencils::BFInputWStencil::applyRightWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputWStencil::applyBottomWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}
void Stencils::BFInputWStencil::applyTopWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j) {}

void Stencils::BFInputWStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  flowField.getTurbOmega().getScalar(i, j, k) = computeW3D(i, j, k, stepSize_, parameters_);
}

void Stencils::BFInputWStencil::applyRightWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputWStencil::applyBottomWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputWStencil::applyTopWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputWStencil::applyFrontWall(
  [[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputWStencil::applyBackWall([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
