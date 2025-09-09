#include "StdAfx.hpp"

#include "GlobalBoundaryFactory.hpp"

GlobalBoundaryFactory::GlobalBoundaryFactory(Parameters& parameters):
  parameters_(parameters) {
  // The parameters will be modified, and therefore are not declared as constants.

  // All stencils are created, disregarding whether they will be used or not. This is less
  // complicated and doesn't seem that costly.

  periodic_[0] = new Stencils::PeriodicBoundaryVelocityStencil(parameters);
  periodic_[1] = new Stencils::PeriodicBoundaryFGHStencil(parameters);
  // periodic_[2] = new Stencils::PeriodicBoundaryTVStencil(parameters);

  moving_[0] = new Stencils::MovingWallVelocityStencil(parameters);
  moving_[1] = new Stencils::MovingWallFGHStencil(parameters);
  // moving_[2] = new Stencils::MovingWallTVStencil(parameters);

  outflow_[0] = new Stencils::NeumannVelocityBoundaryStencil(parameters);
  outflow_[1] = new Stencils::NeumannFGHBoundaryStencil(parameters);
  // outflow_[2] = new Stencils::NeumannTVBoundaryStencil(parameters);

  channelInput_[0] = new Stencils::BFInputVelocityStencil(parameters);
  channelInput_[1] = new Stencils::BFInputFGHStencil(parameters);
  // channelInput_[2] = new Stencils::BFInputTVStencil(parameters);

  // Then, assign them according to the scenario.
  std::string scenario = parameters.simulation.scenario;

  if (scenario == "cavity") {
    // Here, all is about setting the velocity at the boundaries.
    for (int i = 0; i < 6; i++) {
      velocityStencils_[i] = moving_[0];
      FGHStencils_[i]      = moving_[1];
      TVStencils_[i]       = new Stencils::MovingWallTVStencil(parameters); // WS2
      KStencils_[i]        = new Stencils::MovingWallKStencil(parameters);  // WS3
      WStencils_[i]        = new Stencils::MovingWallWStencil(parameters);  // WS3
    }
    parameters.walls.typeLeft   = DIRICHLET;
    parameters.walls.typeRight  = DIRICHLET;
    parameters.walls.typeBottom = DIRICHLET;
    parameters.walls.typeTop    = DIRICHLET;
    parameters.walls.typeFront  = DIRICHLET;
    parameters.walls.typeBack   = DIRICHLET;
  } else if (scenario == "channel") {
    // To the left, we have the input
    velocityStencils_[0] = channelInput_[0];
    FGHStencils_[0]      = channelInput_[1];
    TVStencils_[0]       = new Stencils::BFInputTVStencil(parameters); // WS2
    KStencils_[0]        = new Stencils::BFInputKStencil(parameters);  // WS3
    WStencils_[0]        = new Stencils::BFInputWStencil(parameters);  // WS3


    // To the right, there is an outflow boundary
    velocityStencils_[1] = outflow_[0];
    FGHStencils_[1]      = outflow_[1];
    TVStencils_[1]       = new Stencils::NeumannTVBoundaryStencil(parameters); // WS2
    KStencils_[1]        = new Stencils::NeumannKBoundaryStencil(parameters);  // WS3
    WStencils_[1]        = new Stencils::NeumannWBoundaryStencil(parameters);  // WS3

    // The other walls are moving walls
    for (int i = 2; i < 6; i++) {
      velocityStencils_[i] = moving_[0];
      FGHStencils_[i]      = moving_[1];
      TVStencils_[i]       = new Stencils::MovingWallTVStencil(parameters); // WS2
      KStencils_[i]        = new Stencils::MovingWallKStencil(parameters);  // WS3
      WStencils_[i]        = new Stencils::MovingWallWStencil(parameters);  // WS3
    }
    parameters.walls.typeLeft   = DIRICHLET;
    parameters.walls.typeRight  = NEUMANN;
    parameters.walls.typeBottom = DIRICHLET;
    parameters.walls.typeTop    = DIRICHLET;
    parameters.walls.typeFront  = DIRICHLET;
    parameters.walls.typeBack   = DIRICHLET;
  } else if (scenario == "pressure-channel") {
    // We have Dirichlet conditions for pressure on both sides,
    // hence outflow conditions for the velocities.
    velocityStencils_[0] = outflow_[0];
    FGHStencils_[0]      = outflow_[1];
    TVStencils_[0]       = new Stencils::NeumannTVBoundaryStencil(parameters); // WS2

    // To the right, there is an outflow boundary
    velocityStencils_[1] = outflow_[0];
    FGHStencils_[1]      = outflow_[1];
    TVStencils_[2]       = new Stencils::NeumannTVBoundaryStencil(parameters); // WS2

    // The other walls are moving walls
    for (int i = 2; i < 6; i++) {
      velocityStencils_[i] = moving_[0];
      FGHStencils_[i]      = moving_[1];
      TVStencils_[i]       = new Stencils::MovingWallTVStencil(parameters); // WS2
    }
    parameters.walls.typeLeft   = NEUMANN;
    parameters.walls.typeRight  = NEUMANN;
    parameters.walls.typeBottom = DIRICHLET;
    parameters.walls.typeTop    = DIRICHLET;
    parameters.walls.typeFront  = DIRICHLET;
    parameters.walls.typeBack   = DIRICHLET;
  } else if ((scenario == "periodic-box") || (scenario == "taylor-green")) {
    for (int i = 0; i < 6; i++) {
      velocityStencils_[i] = periodic_[0];
      FGHStencils_[i]      = periodic_[1];
      TVStencils_[i]       = new Stencils::PeriodicBoundaryTVStencil(parameters); // WS2
    }
    parameters.walls.typeLeft   = PERIODIC;
    parameters.walls.typeRight  = PERIODIC;
    parameters.walls.typeBottom = PERIODIC;
    parameters.walls.typeTop    = PERIODIC;
    parameters.walls.typeFront  = PERIODIC;
    parameters.walls.typeBack   = PERIODIC;
  } else {
    throw std::runtime_error("Scenario not recognized");
  }
}

GlobalBoundaryFactory::~GlobalBoundaryFactory() {
  delete moving_[0];
  delete moving_[1];

  delete periodic_[0];
  delete periodic_[1];

  delete outflow_[0];
  delete outflow_[1];

  delete channelInput_[0];
  delete channelInput_[1];

  delete TVStencils_[0];
  delete TVStencils_[1];
  delete TVStencils_[2];
  delete TVStencils_[3];
  delete TVStencils_[4];
  delete TVStencils_[5];

  delete KStencils_[0];
  delete KStencils_[1];
  delete KStencils_[2];
  delete KStencils_[3];
  delete KStencils_[4];
  delete KStencils_[5];

  delete WStencils_[0];
  delete WStencils_[1];
  delete WStencils_[2];
  delete WStencils_[3];
  delete WStencils_[4];
  delete WStencils_[5];
}

GlobalBoundaryIterator<FlowField> GlobalBoundaryFactory::getGlobalBoundaryFGHIterator(FlowField& flowField) {
  if (parameters_.geometry.dim == 2) {
    return GlobalBoundaryIterator<FlowField>(flowField, parameters_, *(FGHStencils_[0]), *(FGHStencils_[1]), *(FGHStencils_[2]), *(FGHStencils_[3]), 1, 0);
  }
  return GlobalBoundaryIterator<FlowField>(
    flowField, parameters_, *(FGHStencils_[0]), *(FGHStencils_[1]), *(FGHStencils_[2]), *(FGHStencils_[3]), *(FGHStencils_[4]), *(FGHStencils_[5]), 1, 0
  );
}

GlobalBoundaryIterator<FlowField> GlobalBoundaryFactory::getGlobalBoundaryVelocityIterator(FlowField& flowField) {
  if (parameters_.geometry.dim == 2) {
    return GlobalBoundaryIterator<FlowField>(flowField, parameters_, *(velocityStencils_[0]), *(velocityStencils_[1]), *(velocityStencils_[2]), *(velocityStencils_[3]), 1, 0);
  }
  return GlobalBoundaryIterator<FlowField>(
    flowField, parameters_, *(velocityStencils_[0]), *(velocityStencils_[1]), *(velocityStencils_[2]), *(velocityStencils_[3]), *(velocityStencils_[4]), *(velocityStencils_[5]), 1, 0
  );
}

GlobalBoundaryIterator<TurbulentFlowField> GlobalBoundaryFactory::getGlobalBoundaryTVIterator(TurbulentFlowField& flowField) {
  if (parameters_.geometry.dim == 2) {
    return GlobalBoundaryIterator<TurbulentFlowField>(flowField, parameters_, *(TVStencils_[0]), *(TVStencils_[1]), *(TVStencils_[2]), *(TVStencils_[3]), 1, 0);
  }
  return GlobalBoundaryIterator<TurbulentFlowField>(
    flowField, parameters_, *(TVStencils_[0]), *(TVStencils_[1]), *(TVStencils_[2]), *(TVStencils_[3]), *(TVStencils_[4]), *(TVStencils_[5]), 1, 0
  );
}

GlobalBoundaryIterator<TurbulentFlowField> GlobalBoundaryFactory::getGlobalBoundaryKIterator(TurbulentFlowField& flowField) {
  if (parameters_.geometry.dim == 2) {
    return GlobalBoundaryIterator<TurbulentFlowField>(flowField, parameters_, *(KStencils_[0]), *(KStencils_[1]), *(KStencils_[2]), *(KStencils_[3]), 1, 0);
  }
  return GlobalBoundaryIterator<TurbulentFlowField>(
    flowField, parameters_, *(KStencils_[0]), *(KStencils_[1]), *(KStencils_[2]), *(KStencils_[3]), *(KStencils_[4]), *(KStencils_[5]), 1, 0
  );
}

GlobalBoundaryIterator<TurbulentFlowField> GlobalBoundaryFactory::getGlobalBoundaryWIterator(TurbulentFlowField& flowField) {
  if (parameters_.geometry.dim == 2) {
    return GlobalBoundaryIterator<TurbulentFlowField>(flowField, parameters_, *(WStencils_[0]), *(WStencils_[1]), *(WStencils_[2]), *(WStencils_[3]), 1, 0);
  }
  return GlobalBoundaryIterator<TurbulentFlowField>(
    flowField, parameters_, *(WStencils_[0]), *(WStencils_[1]), *(WStencils_[2]), *(WStencils_[3]), *(WStencils_[4]), *(WStencils_[5]), 1, 0
  );
}
