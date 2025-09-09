#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"

#include "Solvers/PetscSolver.hpp"
#include "Solvers/SORSolver.hpp"

// TODO most functions are copied from Simulation.cpp. needed changes

// UNCHANGED from Simulation.cpp
TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& turbFlowField):
  Simulation(parameters, turbFlowField),
  turbFlowField_(turbFlowField),
  transportUpdateStencil_(parameters),
  transportUpdateIterator_(turbFlowField_, parameters, transportUpdateStencil_, 1, 0),
  transportRHSStencil_(parameters),
  transportRHSIterator_(turbFlowField_, parameters, transportRHSStencil_, 1, 0),
  maxTVStencil_(parameters),
  maxTVFieldIterator_(turbFlowField, parameters, maxTVStencil_),
  maxTVBoundaryIterator_(turbFlowField, parameters, maxTVStencil_),
  maxKVStencil_(parameters),
  maxKVFieldIterator_(turbFlowField, parameters, maxKVStencil_),
  maxKVBoundaryIterator_(turbFlowField, parameters, maxKVStencil_),
  maxWVStencil_(parameters),
  maxWVFieldIterator_(turbFlowField, parameters, maxWVStencil_),
  maxWVBoundaryIterator_(turbFlowField, parameters, maxWVStencil_),
  tvStencil_(parameters),
  tvIterator_(turbFlowField_, parameters, tvStencil_, 1, 0),
  fghStencil_(parameters),
  fghIterator_(turbFlowField_, parameters, fghStencil_),
  turbObstacleStencil_(parameters),
  turbObstacleIterator_(turbFlowField_, parameters, turbObstacleStencil_),
  globalBoundaryFactory_(parameters),
  wallTVIterator_(globalBoundaryFactory_.getGlobalBoundaryTVIterator(turbFlowField_)),
  wallKIterator_(globalBoundaryFactory_.getGlobalBoundaryKIterator(turbFlowField_)),
  wallWIterator_(globalBoundaryFactory_.getGlobalBoundaryWIterator(turbFlowField_)),
  petscParallelManager_(parameters, turbFlowField_) {}

void TurbulentSimulation::WriteFlowField(int n) {
  // Write the flow field to a file
  std::string   str;
  std::string   filename = "output.log/step" + std::to_string(n) + "rank" + std::to_string(parameters_.parallel.rank) + str + ".txt";
  std::string   line_opening;
  std::ofstream myfile;
  myfile.open(filename);
  int Nx = parameters_.parallel.localSize[0] + 3;
  int Ny = parameters_.parallel.localSize[1] + 3;
  myfile << "Nx: " << Nx << " Ny: " << Ny << std::endl;
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      // if (j == 0 || j == 1 || j == parameters_.parallel.localSize[1] + 2) {
      //   line_opening = "B  ";
      // } else {
      //   line_opening = "I  ";
      // }
      // if (i == 0) {
      //   myfile << line_opening << flowField_.getVelocity().getVector(i, j)[0] << " ";
      // } else {
      //   myfile << flowField_.getVelocity().getVector(i, j)[0] << " ";
      // }
      myfile << flowField_.getPressure().getScalar(i, j) << " ";
    }
    myfile << std::endl;
  }
  // myfile << "   ";
  // for (int k = 0; k < parameters_.parallel.localSize[0] + 3; k++) {
  //   if (k == 0 || k == 1 || k == parameters_.parallel.localSize[0] + 2) {
  //     myfile << "B  ";
  //   } else {
  //     myfile << "I  ";
  //   }
  // }
  myfile << std::endl;
}


void TurbulentSimulation::printFlowField() {
  int Nx = turbFlowField_.getCellsX();
  int Ny = turbFlowField_.getCellsY();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  std::ofstream myfile;
  std::string   filename = "step" + std::to_string(count) + "_rank" + std::to_string(rank) + ".txt";
  // std::string   output   = "RANK" + std::to_string(rank) + "(" + std::to_string(Nx) + ", " + std::to_string(Ny) + ")\n";

  myfile.open(filename);
  int decimals = 10;
  int width    = 20;
  myfile << std::setprecision(decimals);
  myfile << std::fixed;
  myfile.setf(std::ios::left);
  myfile << "RANK" << rank << "  (" << Nx << ", " << Ny << ")\n";
  myfile << "time = " << time << "\n";
  myfile << "--------------NOW PRINTING K-----------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getTurbKineticEnergy().getScalar(i, j) << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING W-----------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getTurbOmega().getScalar(i, j) << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING K RHS-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getTurbKineticEnergyRHS().getScalar(i, j) << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING W RHS-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getTurbOmegaRHS().getScalar(i, j) << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING TV-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getTurbViscosity().getScalar(i, j) << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING u-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getVelocity().getVector(i, j)[0] << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING v-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getVelocity().getVector(i, j)[1] << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING F-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getFGH().getVector(i, j)[0] << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING G-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getFGH().getVector(i, j)[1] << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING RHS-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getRHS().getScalar(i, j) << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------NOW PRINTING P-------------------------\n\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      myfile << std::setw(width) << turbFlowField_.getPressure().getScalar(i, j) << " ";
    }
    myfile << "\n";
  }
  myfile << "\n--------------END OF TIME STEP-------------------------\n";
}

void TurbulentSimulation::printFlowField3D() {
  int Nx = turbFlowField_.getCellsX();
  int Ny = turbFlowField_.getCellsY();
  int Nz = turbFlowField_.getCellsZ();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  std::ofstream myfile;
  std::string   filename = "step" + std::to_string(count) + "_rank" + std::to_string(rank) + ".txt";
  // std::string   output   = "RANK" + std::to_string(rank) + "(" + std::to_string(Nx) + ", " + std::to_string(Ny) + ")\n";

  myfile.open(filename);
  int decimals = 10;
  int width    = 20;
  myfile << std::setprecision(decimals);
  myfile << std::fixed;
  myfile.setf(std::ios::left);
  myfile << "RANK" << rank << "  (" << Nx << ", " << Ny << ")\n";
  myfile << "time = " << time << "\n";
  myfile << "--------------NOW PRINTING K-----------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getTurbKineticEnergy().getScalar(i, j, k) << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING W-----------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getTurbOmega().getScalar(i, j, k) << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING K RHS-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getTurbKineticEnergyRHS().getScalar(i, j, k) << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING W RHS-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getTurbOmegaRHS().getScalar(i, j, k) << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING TV-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getTurbViscosity().getScalar(i, j, k) << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING u-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getVelocity().getVector(i, j, k)[0] << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING v-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getVelocity().getVector(i, j, k)[1] << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING w-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getVelocity().getVector(i, j, k)[2] << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING F-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getFGH().getVector(i, j, k)[0] << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING G-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getFGH().getVector(i, j, k)[1] << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING H-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getFGH().getVector(i, j, k)[2] << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING RHS-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getRHS().getScalar(i, j, k) << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------NOW PRINTING P-------------------------\n\n";
  for (int k = 0; k < Nz; k++) {
    myfile << "k = " << k << "\n";
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << std::setw(width) << turbFlowField_.getPressure().getScalar(i, j, k) << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n--------------END OF TIME STEP-------------------------\n";
}


void TurbulentSimulation::initializeFlowField() {
  if (parameters_.simulation.scenario == "taylor-green") {
    // Currently, a particular initialisation is only required for the taylor-green vortex.
    Stencils::InitTaylorGreenFlowFieldStencil stencil(parameters_);
    FieldIterator<FlowField>                  iterator(turbFlowField_, parameters_, stencil);
    iterator.iterate();
  } else if (parameters_.simulation.scenario == "channel") {
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(turbFlowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
    wallVelocityIterator_.iterateLeft();
    wallVelocityIterator_.iterateRight();
    wallVelocityIterator_.iterateTop();
    wallVelocityIterator_.iterateBottom();

  } else if (parameters_.simulation.scenario == "pressure-channel") {
    // Set pressure boundaries here for left wall
    const RealType value = parameters_.walls.scalarLeft;
    ScalarField&   rhs   = turbFlowField_.getRHS();

    if (parameters_.geometry.dim == 2) {
      const int sizey = turbFlowField_.getNy();
      for (int i = 0; i < sizey + 3; i++) {
        rhs.getScalar(0, i) = value;
      }
    } else {
      const int sizey = turbFlowField_.getNy();
      const int sizez = turbFlowField_.getNz();
      for (int i = 0; i < sizey + 3; i++) {
        for (int j = 0; j < sizez + 3; j++) {
          rhs.getScalar(0, i, j) = value;
        }
      }
    }

    // Do same procedure for domain flagging as for regular channel
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(turbFlowField_, parameters_, stencil);
    iterator.iterate();
  }

  // WS2 wall distance iterator
  Stencils::WallDistanceStencil     stencil(parameters_);
  FieldIterator<TurbulentFlowField> iterator(turbFlowField_, parameters_, stencil, 1, 0);
  iterator.iterate();

  // WS3 initialize k and omega
  Stencils::InitKWStencil           IKWstencil(parameters_);
  FieldIterator<TurbulentFlowField> IKWiterator(turbFlowField_, parameters_, IKWstencil, 1, 0);
  IKWiterator.iterate();
  wallKIterator_.iterateTop();
  wallKIterator_.iterateBottom();
  wallKIterator_.iterateRight();
  wallKIterator_.iterateLeft();

  wallWIterator_.iterateTop();
  wallWIterator_.iterateBottom();
  wallWIterator_.iterateRight();
  wallWIterator_.iterateLeft();

  wallTVIterator_.iterateTop();
  wallTVIterator_.iterateBottom();
  wallTVIterator_.iterateRight();
  wallTVIterator_.iterateLeft();


  solver_->reInitMatrix();
}

void TurbulentSimulation::solveTimestep() {
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();

  // printFlowField3D();
  // count++;

  transportRHSIterator_.iterate();    // compute kw model rhs terms
  transportUpdateIterator_.iterate(); // compute kw values

  wallKIterator_.iterateTop();
  wallKIterator_.iterateBottom();
  wallKIterator_.iterateRight();
  wallKIterator_.iterateLeft(); // k&w global boundary values

  wallWIterator_.iterateTop();
  wallWIterator_.iterateBottom();
  wallWIterator_.iterateRight();
  wallWIterator_.iterateLeft();

  tvIterator_.iterate(); // compute tv from kw
  wallTVIterator_.iterateTop();
  wallTVIterator_.iterateBottom();
  wallTVIterator_.iterateRight();
  wallTVIterator_.iterateLeft(); // tv global boundary values, computed from kw


  petscParallelManager_.communicateTurbViscosity();
  petscParallelManager_.communicateTurbKineticEnergy();
  petscParallelManager_.communicateTurbOmega();
  //  Compute FGH
  fghIterator_.iterate();
  // Set global boundary values
  wallFGHIterator_.iterateTop();
  wallFGHIterator_.iterateBottom();
  wallFGHIterator_.iterateLeft();
  wallFGHIterator_.iterateRight();
  // compute the right hand side(RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();
  // Communicate pressure
  petscParallelManager_.communicatePressure();
  // Compute velocity
  velocityIterator_.iterate();
  turbObstacleIterator_.iterate();
  // Communicate velocity
  petscParallelManager_.communicateVelocity();

  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterateLeft();
  wallVelocityIterator_.iterateRight();
  wallVelocityIterator_.iterateTop();
  wallVelocityIterator_.iterateBottom();

  // printFlowField3D();
  // std::cin.get();
  // count++;
  // time += parameters_.timestep.dt;
}

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::TurbVTKStencil          turbVtkStencil(parameters_);
  FieldIterator<TurbulentFlowField> turbVtkIterator(turbFlowField_, parameters_, turbVtkStencil, 1, 0);

  turbVtkIterator.iterate();
  turbVtkStencil.write(turbFlowField_, timeStep, simulationTime);
}

void TurbulentSimulation::setTimeStep() {
  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin()) + 1.0 / (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin());
  // Determine maximum velocity
  maxTVStencil_.reset();
  maxTVFieldIterator_.iterate();
  maxTVBoundaryIterator_.iterate();
  maxUStencil_.reset();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();
  if (parameters_.geometry.dim == 3) {
    factor += 1.0 / (parameters_.meshsize->getDzMin() * parameters_.meshsize->getDzMin());
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2] + std::numeric_limits<RealType>::epsilon());
  } else {
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0] + std::numeric_limits<RealType>::epsilon());
  }
  // std::cout << *maxTVStencil_.getMaxValues() << std::endl;
  // introduce the turbolent viscosity in the timestep calculation)

  // localMin = std::min(
  //   std::min(parameters_.flow.Re / (2 * factor), (1.0 / (maxTVStencil_.getMaxValues()[0] + 1.0 / parameters_.flow.Re)) / (2 * factor)),
  //   std::min(parameters_.timestep.dt, std::min(1 / (maxUStencil_.getMaxValues()[0] + MY_EPS), 1 / (maxUStencil_.getMaxValues()[1] + MY_EPS)))
  // )

  localMin = std::min(
    std::min(parameters_.flow.Re / (2 * factor), (1.0 / (maxTVStencil_.getMaxValues()[0] + 1.0 / parameters_.flow.Re)) / (2 * factor)),
    std::min(
      std::min(parameters_.timestep.dt, std::min(1 / (maxUStencil_.getMaxValues()[0] + MY_EPS), 1 / (maxUStencil_.getMaxValues()[1] + MY_EPS))),
      std::min((1.0 / (maxWVStencil_.getMaxValues()[0] + MY_EPS)) / (2 * factor), (1.0 / (maxKVStencil_.getMaxValues()[0] + MY_EPS)) / (2 * factor))
    )
  );


  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}
