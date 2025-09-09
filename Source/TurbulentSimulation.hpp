#pragma once

#include "ParallelManagers/PetscParallelManagerTurbulent.hpp"
#include "Simulation.hpp"
#include "Stencils/InitKWStencil.hpp"
#include "Stencils/MaxKVStencil.hpp"
#include "Stencils/MaxTVStencil.hpp"
#include "Stencils/MaxWVStencil.hpp"
#include "Stencils/TransportRHSStencil.hpp"
#include "Stencils/TransportUpdateStencil.hpp"
#include "Stencils/TurbVTKStencil.hpp"
#include "Stencils/TVStencil.hpp"
#include "Stencils/WallDistanceStencil.hpp"

class TurbulentSimulation: public Simulation { // WS2 TODO 5.5.2
protected:
  TurbulentFlowField& turbFlowField_;

  Stencils::TransportUpdateStencil  transportUpdateStencil_;
  FieldIterator<TurbulentFlowField> transportUpdateIterator_;

  Stencils::TransportRHSStencil     transportRHSStencil_;
  FieldIterator<TurbulentFlowField> transportRHSIterator_;

  Stencils::MaxTVStencil                     maxTVStencil_;
  FieldIterator<TurbulentFlowField>          maxTVFieldIterator_;
  GlobalBoundaryIterator<TurbulentFlowField> maxTVBoundaryIterator_;

  Stencils::MaxKVStencil                     maxKVStencil_;
  FieldIterator<TurbulentFlowField>          maxKVFieldIterator_;
  GlobalBoundaryIterator<TurbulentFlowField> maxKVBoundaryIterator_;

  Stencils::MaxWVStencil                     maxWVStencil_;
  FieldIterator<TurbulentFlowField>          maxWVFieldIterator_;
  GlobalBoundaryIterator<TurbulentFlowField> maxWVBoundaryIterator_;

  Stencils::TVStencil               tvStencil_;
  FieldIterator<TurbulentFlowField> tvIterator_;

  Stencils::TurbFGHStencil          fghStencil_;
  FieldIterator<TurbulentFlowField> fghIterator_;

  Stencils::TurbObstacleStencil     turbObstacleStencil_;
  FieldIterator<TurbulentFlowField> turbObstacleIterator_;

  // Set up the boundary conditions
  GlobalBoundaryFactory                      globalBoundaryFactory_;
  GlobalBoundaryIterator<TurbulentFlowField> wallTVIterator_;
  GlobalBoundaryIterator<TurbulentFlowField> wallKIterator_;
  GlobalBoundaryIterator<TurbulentFlowField> wallWIterator_;

  ParallelManagers::PetscParallelManagerTurbulent petscParallelManager_;

  virtual void setTimeStep() override;

public:
  unsigned int count;
  RealType     time;
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& turbFlowField);
  virtual ~TurbulentSimulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField() override;

  virtual void solveTimestep() override;

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime) override;

  virtual void printFlowField();
  virtual void printFlowField3D();

  virtual void WriteFlowField(int n);
};