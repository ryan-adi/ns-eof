#include "StdAfx.hpp"

#include "InitKWStencil.hpp"

/**
 * InitKWStencil
 *
 */

Stencils::InitKWStencil::InitKWStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  xLimit_(parameters.bfStep.xRatio * parameters.geometry.lengthX),
  yLimit_(parameters.bfStep.yRatio * parameters.geometry.lengthY) {}

void Stencils::InitKWStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const int      obstacle      = flowField.getFlags().getValue(i, j);
  const RealType Re_dh         = parameters_.flow.Re * (1.0 - parameters_.bfStep.yRatio) / parameters_.geometry.lengthY;
  const RealType turbIntensity = 0.16 * std::pow(Re_dh, -0.125);
  const RealType turbLenght    = 0.07 * parameters_.geometry.lengthY;                          // 7% of the hydraulic diameter, for channel flow it is the height of the channel
  const RealType turbKinEnergy = 3.0 / 2.0 * turbIntensity * turbIntensity;                    // 3/2 * (U*I)^2 with U = 1.0
  const RealType turbOmega     = std::pow(0.09, -0.25) * std::sqrt(turbKinEnergy) / turbLenght; // C_mu^1/4 * sqrt(k) / l
  const RealType tv            = turbKinEnergy / (turbOmega + MY_EPS);                         // k / (omega + eps)

  if ((obstacle & OBSTACLE_SELF) == 0) {
    flowField.getTurbKineticEnergy().getScalar(i, j) = turbKinEnergy;
    flowField.getTurbOmega().getScalar(i, j)         = turbOmega;
    flowField.getTurbViscosity().getScalar(i, j)     = tv;
  } else {
    flowField.getTurbKineticEnergy().getScalar(i, j) = 0.0;
    flowField.getTurbOmega().getScalar(i, j)         = 0.0;
    flowField.getTurbViscosity().getScalar(i, j)     = 0.0;
  }
}

// 3d case
void Stencils::InitKWStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  const int      obstacle      = flowField.getFlags().getValue(i, j, k);
  const RealType Re_dh         = parameters_.flow.Re * (1.0 - parameters_.bfStep.yRatio) / parameters_.geometry.lengthY;
  const RealType turbIntensity = 0.16 * std::pow(Re_dh, -0.125);
  const RealType turbLenght    = 0.07 * parameters_.geometry.lengthY;                          // 7% of the hydraulic diameter, for channel flow it is the height of the channel
  const RealType turbKinEnergy = 3.0 / 2.0 * turbIntensity * turbIntensity;                    // 3/2 * (U*I)^2 with U = 1.0
  const RealType turbOmega     = std::pow(0.09, -0.25) * std::sqrt(turbKinEnergy) / turbLenght; // C_mu^1/4 * sqrt(k) / l
  const RealType tv            = turbKinEnergy / (turbOmega + MY_EPS);                         // k / (omega + eps)

  if ((obstacle & OBSTACLE_SELF) == 0) {
    flowField.getTurbKineticEnergy().getScalar(i, j, k) = turbKinEnergy;
    flowField.getTurbOmega().getScalar(i, j, k)         = turbOmega;
    flowField.getTurbViscosity().getScalar(i, j, k)     = tv;
  } else {
    flowField.getTurbKineticEnergy().getScalar(i, j, k) = 0.0;
    flowField.getTurbOmega().getScalar(i, j, k)         = 0.0;
    flowField.getTurbViscosity().getScalar(i, j, k)     = 0.0;
  }
};
