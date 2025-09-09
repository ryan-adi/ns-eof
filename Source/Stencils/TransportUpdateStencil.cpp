#include "StdAfx.hpp"

#include "TransportUpdateStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"


Stencils::TransportUpdateStencil::TransportUpdateStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::TransportUpdateStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  RealType& K      = flowField.getTurbKineticEnergy().getScalar(i, j);
  RealType& RHS_K  = flowField.getTurbKineticEnergyRHS().getScalar(i, j);
  RealType& W      = flowField.getTurbOmega().getScalar(i, j);
  RealType& RHS_W  = flowField.getTurbOmegaRHS().getScalar(i, j);
  RealType& K_prev = flowField.getPrevTurbKineticEnergy().getScalar(i, j);
  RealType& W_prev = flowField.getPrevTurbOmega().getScalar(i, j);
  RealType  dt     = parameters_.timestep.dt;

  // Second order update
  RealType K_current = K;
  RealType W_current = W;

  K = (4 * K - K_prev + 2 * dt * RHS_K) / 3;
  W = (4 * W - W_prev + 2 * dt * RHS_W) / 3;

  K_prev = K_current;
  W_prev = W_current;

  // First order update
  // K = K + parameters_.timestep.dt * RHS_K;
  // W = W + parameters_.timestep.dt * RHS_W;
}

void Stencils::TransportUpdateStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType& K      = flowField.getTurbKineticEnergy().getScalar(i, j, k);
  RealType& RHS_K  = flowField.getTurbKineticEnergyRHS().getScalar(i, j, k);
  RealType& W      = flowField.getTurbOmega().getScalar(i, j, k);
  RealType& RHS_W  = flowField.getTurbOmegaRHS().getScalar(i, j, k);
  RealType& K_prev = flowField.getPrevTurbKineticEnergy().getScalar(i, j, k);
  RealType& W_prev = flowField.getPrevTurbOmega().getScalar(i, j, k);
  RealType  dt     = parameters_.timestep.dt;

  // Second order update
  RealType K_current = K;
  RealType W_current = W;

  K = (4 * K - K_prev + 2 * dt * RHS_K) / 3;
  W = (4 * W - W_prev + 2 * dt * RHS_W) / 3;

  K_prev = K_current;
  W_prev = W_current;

  // First order update
  // K = K + parameters_.timestep.dt * RHS_K;
  // W = W + parameters_.timestep.dt * RHS_W;
}
