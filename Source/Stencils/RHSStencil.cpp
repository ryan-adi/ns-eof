#include "StdAfx.hpp"

#include "RHSStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::RHSStencil::RHSStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j) {
  // Load local velocities into the center layer of the local array
  // loadLocalVelocity2D(flowField, localVelocity_, i, j);
  // loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

  // need to compute dFdx, dGdy
  RealType* const FG       = flowField.getFGH().getVector(i, j);
  RealType* const FG_prevx = flowField.getFGH().getVector(i - 1, j);
  RealType* const FG_prevy = flowField.getFGH().getVector(i, j - 1);

  const RealType dFdx = (FG[0] - FG_prevx[0]) / parameters_.meshsize->getDx(i, j);
  const RealType dGdy = (FG[1] - FG_prevy[1]) / parameters_.meshsize->getDy(i, j);

  flowField.getRHS().getScalar(i, j) = (dFdx + dGdy) / parameters_.timestep.dt;
  // std::cout << "RHSStencil working on cell (" << i << ", " << j << ")\n";
}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j, int k) {
  // The same as in 2D, with slight modifications.

  const int obstacle = flowField.getFlags().getValue(i, j, k);

  if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
    // loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
    // loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);
    RealType* const FGH       = flowField.getFGH().getVector(i, j, k);
    RealType* const FGH_prevx = flowField.getFGH().getVector(i - 1, j, k);
    RealType* const FGH_prevy = flowField.getFGH().getVector(i, j - 1, k);
    RealType* const FGH_prevz = flowField.getFGH().getVector(i, j, k - 1);

    const RealType dFdx = (FGH[0] - FGH_prevx[0]) / parameters_.meshsize->getDx(i, j, k);
    const RealType dGdy = (FGH[1] - FGH_prevy[1]) / parameters_.meshsize->getDy(i, j, k);
    const RealType dHdz = (FGH[2] - FGH_prevz[2]) / parameters_.meshsize->getDz(i, j, k);

    flowField.getRHS().getScalar(i, j, k) = (dFdx + dGdy + dHdz) / parameters_.timestep.dt;
    // if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
    // }
    // if ((obstacle & OBSTACLE_TOP) == 0) {
    //   values[1] = computeG3D(localVelocity_, localMeshsize_, parameters_, parameters_.timestep.dt);
    // }
    // if ((obstacle & OBSTACLE_BACK) == 0) {
    //   values[2] = computeH3D(localVelocity_, localMeshsize_, parameters_, parameters_.timestep.dt);
    // }
  }
}
