#include "StdAfx.hpp"

#include "TVStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

/**
 * TVStencil apply 2D mostly done needs recheck, 3D analog to 2D
 *
 */

Stencils::TVStencil::TVStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::TVStencil::apply(TurbulentFlowField& flowField, int i, int j) {

  if (parameters_.simulation.type == "turbulence_ml") {
    // VectorField&   velocity = flowField.getVelocity();
    //  Load local velocities into the center layer of the local array RA: TODO use local velocity ???
    loadLocalVelocity2D(flowField, localVelocity_, i, j);
    loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);
    std::string    bl_type       = parameters_.turb.bl_type; // should only be "non-existing", "lam-flat", "turb-flat"
    const RealType Re            = parameters_.flow.Re;
    const RealType posX          = parameters_.meshsize->getPosX(i, j); // x-pos left bottom front node
    const RealType dx            = parameters_.meshsize->getDx(i, j);
    ScalarField&   nearestWallH  = flowField.getNearestWallH();
    ScalarField&   turbViscosity = flowField.getTurbViscosity();

    // calculate Re_x
    const RealType Re_x = Re * (posX + 0.5 * dx) + MY_EPS; //* velocity.getVector(i, j)[0] + MY_EPS;

    // define delta based on boundary layer type
    RealType delta = 0;

    if (bl_type == "non-existing") {
      delta = 0;
    } else if (bl_type == "lam-flat") {
      delta = 4.91 * (posX + 0.5 * dx) / std::pow(Re_x, 0.5);
    } else if (bl_type == "turb-flat") {
      delta = 0.382 * (posX + 0.5 * dx) / std::pow(Re_x, 0.2);
    } else { // TODO move error to xml read command
      throw std::runtime_error("Boundary layer type not recognized");
    };

    // calculate mixing length
    const RealType l_m = std::min(0.41 * nearestWallH.getScalar(i, j), 0.09 * delta);

    // calculate turbulent viscosity
    const RealType SQuad          = computeSQuad2D(localVelocity_, localMeshsize_);
    turbViscosity.getScalar(i, j) = std::pow(l_m, 2) * std::pow(2 * SQuad, 0.5);
  } else if (parameters_.simulation.type == "turbulence_kw") {
    RealType& turbViscosity = flowField.getTurbViscosity().getScalar(i, j);

    turbViscosity = flowField.getTurbKineticEnergy().getScalar(i, j) / (flowField.getTurbOmega().getScalar(i, j) + MY_EPS);
  } else {
    throw std::runtime_error("Unknown simulation type! Currently supported: turbulence_ml, turbulence_kw");
  }
}

void Stencils::TVStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {

  if (parameters_.simulation.type == "turbulence_ml") {
    // VectorField&   velocity = flowField.getVelocity();
    //  Load local velocities into the center layer of the local array RA: TODO use local velocity ???
    loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
    loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);
    std::string    bl_type       = parameters_.turb.bl_type; // should only be "non-existing", "lam-flat", "turb-flat"
    const RealType Re            = parameters_.flow.Re;
    const RealType posX          = parameters_.meshsize->getPosX(i, j, k);
    const RealType dx            = parameters_.meshsize->getDx(i, j, k);
    ScalarField&   nearestWallH  = flowField.getNearestWallH();
    ScalarField&   turbViscosity = flowField.getTurbViscosity();

    // calculate Re_x
    const RealType Re_x = Re * (posX + 0.5 * dx) + MY_EPS;

    // define delta based on boundary layer type
    RealType delta = 0;
    if (bl_type == "non-existing") {
      delta = 0;
    } else if (bl_type == "lam-flat") {
      delta = 4.91 * (posX + 0.5 * dx) / std::pow(Re_x, 0.5);
    } else if (bl_type == "turb-flat") {
      delta = 0.382 * (posX + 0.5 * dx) / std::pow(Re_x, 0.2);
    } else {
      throw std::runtime_error("Boundary layer type not recognized");
    };

    // calculate mixing lnegth
    const RealType l_m = std::min(0.41 * nearestWallH.getScalar(i, j, k), 0.09 * delta);

    // calculate turbulent viscosity
    const RealType SQuad             = computeSQuad3D(localVelocity_, localMeshsize_);
    turbViscosity.getScalar(i, j, k) = std::pow(l_m, 2) * std::pow(2 * SQuad, 0.5);
  } else if (parameters_.simulation.type == "turbulence_kw") {

    RealType& turbViscosity = flowField.getTurbViscosity().getScalar(i, j, k);
    turbViscosity           = flowField.getTurbKineticEnergy().getScalar(i, j, k) / flowField.getTurbOmega().getScalar(i, j, k);
  } else {
    throw std::runtime_error("Unknown simulation type! Currently supported: turbulence_ml, turbulence_kw");
  }
}
