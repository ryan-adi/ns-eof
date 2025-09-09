#include "StdAfx.hpp"

#include "Simulation.hpp"

#include <fstream>
#include <iostream>

#include "Solvers/PetscSolver.hpp"
#include "Solvers/SORSolver.hpp"

Simulation::Simulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField),
  maxUStencil_(parameters),
  maxUFieldIterator_(flowField_, parameters, maxUStencil_),
  maxUBoundaryIterator_(flowField_, parameters, maxUStencil_),
  globalBoundaryFactory_(parameters),
  wallVelocityIterator_(globalBoundaryFactory_.getGlobalBoundaryVelocityIterator(flowField_)),
  wallFGHIterator_(globalBoundaryFactory_.getGlobalBoundaryFGHIterator(flowField_)),
  fghStencil_(parameters),
  fghIterator_(flowField_, parameters, fghStencil_),
  rhsStencil_(parameters),
  rhsIterator_(flowField_, parameters, rhsStencil_),
  velocityStencil_(parameters),
  obstacleStencil_(parameters),
  velocityIterator_(flowField_, parameters, velocityStencil_),
  obstacleIterator_(flowField_, parameters, obstacleStencil_),
  petscParallelManager_(parameters, flowField_)
#ifdef ENABLE_PETSC
  ,
  solver_(std::make_unique<Solvers::PetscSolver>(flowField_, parameters))
#else
  ,
  solver_(std::make_unique<Solvers::SORSolver>(flowField_, parameters))
#endif
{
}

void Simulation::printFlowField() {
  int Nx = flowField_.getCellsX();
  int Ny = flowField_.getCellsY();

  // int rank;
  // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // std::string flag_str = flag == 0 ? "Before com:" : "After com:";

  std::string output = "RANK ? (" + std::to_string(Nx) + ", " + std::to_string(Ny) + ")\n";

  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      std::string val = std::to_string(flowField_.getVelocity().getVector(i, j)[0]);
      // std::string val = std::to_string(flowField_.getFlags().getValue(i,j));
      output += val + " ";
    }
    output += "\n";
  }
  output += "----------------------------------------------------------\n";
  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      std::string val = std::to_string(flowField_.getFGH().getVector(i, j)[0]);
      // std::string val = std::to_string(turbFlowField_.getFlags().getValue(i,j));
      output += val + " ";
    }
    output += "\n";
  }

  std::cout << "-----------------------------------\n" << output << "-----------------------------------\n";
}

void Simulation::initializeFlowField() {
  if (parameters_.simulation.scenario == "taylor-green") {
    // Currently, a particular initialisation is only required for the taylor-green vortex.
    Stencils::InitTaylorGreenFlowFieldStencil stencil(parameters_);
    FieldIterator<FlowField>                  iterator(flowField_, parameters_, stencil);
    iterator.iterate();
  } else if (parameters_.simulation.scenario == "channel") {
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(flowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
    wallVelocityIterator_.iterate();
  } else if (parameters_.simulation.scenario == "pressure-channel") {
    // Set pressure boundaries here for left wall
    const RealType value = parameters_.walls.scalarLeft;
    ScalarField&   rhs   = flowField_.getRHS();

    if (parameters_.geometry.dim == 2) {
      const int sizey = flowField_.getNy();
      for (int i = 0; i < sizey + 3; i++) {
        rhs.getScalar(0, i) = value;
      }
    } else {
      const int sizey = flowField_.getNy();
      const int sizez = flowField_.getNz();
      for (int i = 0; i < sizey + 3; i++) {
        for (int j = 0; j < sizez + 3; j++) {
          rhs.getScalar(0, i, j) = value;
        }
      }
    }

    // Do same procedure for domain flagging as for regular channel
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(flowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
  }

  solver_->reInitMatrix();
}
void Simulation::fillFlowField() {
  int Nx = flowField_.getCellsX();
  int Ny = flowField_.getCellsY();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++) {
      // flowField_.getVelocity().getVector(i, j)[0] = rank;
      // flowField_.getVelocity().getVector(i, j)[1] = -1 * rank;
      flowField_.getPressure().getScalar(i, j) = rank;
    }
  }
}


void Simulation::fillFlowField3D() {
  int Nx = flowField_.getCellsX();
  int Ny = flowField_.getCellsY();
  int Nz = flowField_.getCellsZ();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (int j = Ny - 1; j >= 0; j--) {
    for (int i = 0; i < Nx; i++)
      for (int k = 0; k < Nz; k++) {
        // flowField_.getVelocity().getVector(i, j)[0] = rank;
        // flowField_.getVelocity().getVector(i, j)[1] = -1 * rank;
        flowField_.getPressure().getScalar(i, j, k) = rank;
      }
  }
}
void Simulation::WriteFlowField(int n, int flag) {
  // Write the flow field to a file
  std::string comm;
  if (flag == 0) {
    comm = "BC";
  } else if (flag == 1) {
    comm = "AC";
  }
  std::string var_name;
  for (int k = 0; k < 1; k++) {

    // if (k == 0) {
    //   var_name = "u";
    // } else if (k == 1) {
    //   var_name = "v";
    // }
    std::string   filename = "step_" + std::to_string(n) + "_" + "p" + "_rank" + std::to_string(parameters_.parallel.rank) + "_" + comm + ".txt";
    std::ofstream myfile;
    myfile.open(filename);
    int Nx = parameters_.parallel.localSize[0] + 3;
    int Ny = parameters_.parallel.localSize[1] + 3;
    myfile << "Nx: " << Nx << " Ny: " << Ny << std::endl;
    for (int j = Ny - 1; j >= 0; j--) {
      for (int i = 0; i < Nx; i++) {
        myfile << flowField_.getPressure().getScalar(i, j) << " ";
      }
      myfile << std::endl;
    }
    myfile << std::endl;
    myfile.close();
  }
}

void Simulation::WriteFlowField3D(int n, int flag) {
  // Write the flow field to a file
  std::string comm;
  if (flag == 0) {
    comm = "BC";
  } else if (flag == 1) {
    comm = "AC";
  }
  std::string var_name;
  for (int kk = 0; kk < 1; kk++) {

    // if (k == 0) {
    //   var_name = "u";
    // } else if (k == 1) {
    //   var_name = "v";
    // }
    std::string   filename = "output.log/step_" + std::to_string(n) + "_" + "p" + "_rank" + std::to_string(parameters_.parallel.rank) + "_" + comm + ".txt";
    std::ofstream myfile;
    myfile.open(filename);
    int Nx = parameters_.parallel.localSize[0] + 3;
    int Ny = parameters_.parallel.localSize[1] + 3;
    int Nz = parameters_.parallel.localSize[2] + 3;
    myfile << "Nx: " << Nx << " Ny: " << Ny << std::endl;
    for (int k = 0; k < Nz; k++) {
      for (int j = Ny - 1; j >= 0; j--) {
        for (int i = 0; i < Nx; i++) {

          myfile << flowField_.getPressure().getScalar(i, j, k) << " ";
        }
        myfile << std::endl;
      }
      myfile << std::endl;
      myfile << std::endl;
      myfile << std::endl;
    }


    myfile.close();
  }
}

void Simulation::solveTimestep() {
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();

  // Compute FGH
  fghIterator_.iterate();
  // Set global boundary values
  // printFlowField();
  wallFGHIterator_.iterate();
  // compute the right hand side(RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();
  // communicate pressure values
  petscParallelManager_.communicatePressure();
  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();
  // communicate velocity values
  petscParallelManager_.communicateVelocity();
  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();
}

void Simulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::VTKStencil     vtkStencil(parameters_);
  FieldIterator<FlowField> vtkIterator(flowField_, parameters_, vtkStencil, 1, 0);

  vtkIterator.iterate();
  vtkStencil.write(flowField_, timeStep, simulationTime);
}

void Simulation::setTimeStep() {
  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin()) + 1.0 / (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin());
  // Determine maximum velocity
  maxUStencil_.reset();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();
  if (parameters_.geometry.dim == 3) {
    factor += 1.0 / (parameters_.meshsize->getDzMin() * parameters_.meshsize->getDzMin());
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2] + std::numeric_limits<RealType>::epsilon());
  } else {
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0] + std::numeric_limits<RealType>::epsilon());
  }

  localMin = std::min(
    parameters_.flow.Re / (2 * factor),
    std::min(
      parameters_.timestep.dt,
      std::min(1 / (maxUStencil_.getMaxValues()[0] + std::numeric_limits<RealType>::epsilon()), 1 / (maxUStencil_.getMaxValues()[1] + std::numeric_limits<RealType>::epsilon()))
    )
  );

  // std::cout
  //   << "x-comp " << maxUStencil_.getMaxValues()[0] << std::endl
  //   << "y-comp " << maxUStencil_.getMaxValues()[1] << std::endl
  //   << "z-comp " << maxUStencil_.getMaxValues()[2] << std::endl;

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}
