#include "StdAfx.hpp"

#include "ObstacleStencil.hpp"

Stencils::ObstacleStencil::ObstacleStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::ObstacleStencil::apply(FlowField& flowField, int i, int j) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();

  // Check if current cell is obstacle cell
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // If top cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_TOP) == 0) {
      const RealType dy_t         = parameters_.meshsize->getDy(i, j + 1);
      const RealType dy           = parameters_.meshsize->getDy(i, j);
      velocity.getVector(i, j)[0] = -dy / dy_t * velocity.getVector(i, j + 1)[0];
    }
    // Same for bottom
    if ((obstacle & OBSTACLE_BOTTOM) == 0) {
      const RealType dy_b         = parameters_.meshsize->getDy(i, j - 1);
      const RealType dy           = parameters_.meshsize->getDy(i, j);
      velocity.getVector(i, j)[0] = -dy / dy_b * velocity.getVector(i, j - 1)[0];
    }
    // If right cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      const RealType dx_r         = parameters_.meshsize->getDx(i + 1, j);
      const RealType dx           = parameters_.meshsize->getDx(i, j);
      velocity.getVector(i, j)[1] = -dx / dx_r * velocity.getVector(i + 1, j)[1];
    }
    // Same for left
    if ((obstacle & OBSTACLE_LEFT) == 0) {
      const RealType dx_l         = parameters_.meshsize->getDx(i - 1, j);
      const RealType dx           = parameters_.meshsize->getDx(i, j);
      velocity.getVector(i, j)[1] = -dx / dx_l * velocity.getVector(i - 1, j)[1];
    }

    // Set normal velocity to zero if right neighbour is not obstacle
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      velocity.getVector(i, j)[0] = 0.0;
    }

    // Set normal velocity to zero if top neighbour is not obstacle
    if ((obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVector(i, j)[1] = 0.0;
    }
  }
}

void Stencils::ObstacleStencil::apply(FlowField& flowField, int i, int j, int k) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();

  // Check if current cell is obstacle cell
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // If top cell is fluid: two velocities have to be set: direction 0 and 2.
    if ((obstacle & OBSTACLE_TOP) == 0) {
      const RealType dy_t            = parameters_.meshsize->getDy(i, j + 1, k);
      const RealType dy              = parameters_.meshsize->getDy(i, j, k);
      velocity.getVector(i, j, k)[0] = -dy / dy_t * velocity.getVector(i, j + 1, k)[0];
      velocity.getVector(i, j, k)[2] = -dy / dy_t * velocity.getVector(i, j + 1, k)[2];
    }
    if ((obstacle & OBSTACLE_BOTTOM) == 0) {
      const RealType dy_b            = parameters_.meshsize->getDy(i, j - 1, k);
      const RealType dy              = parameters_.meshsize->getDy(i, j, k);
      velocity.getVector(i, j, k)[0] = -dy / dy_b * velocity.getVector(i, j - 1, k)[0];
      velocity.getVector(i, j, k)[2] = -dy / dy_b * velocity.getVector(i, j - 1, k)[2];
    }

    // If right cell is fluid: two velocities have to be set: direction 1 and 2.
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      const RealType dx_r            = parameters_.meshsize->getDx(i + 1, j, k);
      const RealType dx              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dx / dx_r * velocity.getVector(i + 1, j, k)[1];
      velocity.getVector(i, j, k)[2] = -dx / dx_r * velocity.getVector(i + 1, j, k)[2];
    }
    if ((obstacle & OBSTACLE_LEFT) == 0) {
      const RealType dx_l            = parameters_.meshsize->getDx(i - 1, j, k);
      const RealType dx              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dx / dx_l * velocity.getVector(i - 1, j, k)[1];
      velocity.getVector(i, j, k)[2] = -dx / dx_l * velocity.getVector(i - 1, j, k)[2];
    }

    // Same for fluid cell in front
    if ((obstacle & OBSTACLE_BACK) == 0) {
      const RealType dz_f            = parameters_.meshsize->getDx(i, j, k + 1);
      const RealType dz              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dz / dz_f * velocity.getVector(i, j, k + 1)[1];
      velocity.getVector(i, j, k)[0] = -dz / dz_f * velocity.getVector(i, j, k + 1)[0];
    }
    if ((obstacle & OBSTACLE_FRONT) == 0) {
      const RealType dz_b            = parameters_.meshsize->getDx(i, j, k - 1);
      const RealType dz              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dz / dz_b * velocity.getVector(i, j, k - 1)[1];
      velocity.getVector(i, j, k)[0] = -dz / dz_b * velocity.getVector(i, j, k - 1)[0];
    }

    // Now the normal velocities need to be set to zero to ensure no flow at interfaces between solid and fluid.
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      velocity.getVector(i, j, k)[0] = 0.0;
    }
    if ((obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVector(i, j, k)[1] = 0.0;
    }
    if ((obstacle & OBSTACLE_BACK) == 0) {
      velocity.getVector(i, j, k)[2] = 0.0;
    }
  }
}

RealType Stencils::TurbObstacleStencil::computeW2D(RealType Re, RealType Dx) {
  RealType beta1 = 3.0 / 40.0;
  return 6.0 / (Re * beta1 * Dx * Dx);
}

Stencils::TurbObstacleStencil::TurbObstacleStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::TurbObstacleStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();
  RealType&    TK       = flowField.getTurbKineticEnergy().getScalar(i, j);
  RealType&    W        = flowField.getTurbOmega().getScalar(i, j);
  RealType     Re       = parameters_.flow.Re;

  // Check if current cell is obstacle cell
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // If top cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_TOP) == 0) {
      const RealType dy_t         = parameters_.meshsize->getDy(i, j + 1);
      const RealType dy           = parameters_.meshsize->getDy(i, j);
      velocity.getVector(i, j)[0] = -dy / dy_t * velocity.getVector(i, j + 1)[0];
      
      RealType omega_w = computeW2D(Re, dy_t / 2);
      RealType omega_1 = flowField.getTurbOmega().getScalar(i, j + 1);
      W                = ((dy + dy_t) / dy_t) * (omega_w - omega_1 * (dy / (dy + dy_t)));
      TK               = -dy / dy_t * flowField.getTurbKineticEnergy().getScalar(i, j + 1);
    }

    // Bottom cell is a fluid cell
    if ((obstacle & OBSTACLE_BOTTOM) == 0) {
      const RealType dy_b         = parameters_.meshsize->getDy(i, j - 1);
      const RealType dy           = parameters_.meshsize->getDy(i, j);
      velocity.getVector(i, j)[0] = -dy / dy_b * velocity.getVector(i, j - 1)[0];
      RealType omega_w = computeW2D(Re, dy_b / 2);
      RealType omega_1 = flowField.getTurbOmega().getScalar(i, j - 1);
      W                = ((dy + dy_b) / dy_b) * (omega_w - omega_1 * (dy / (dy + dy_b)));
      TK                          = -dy / dy_b * flowField.getTurbKineticEnergy().getScalar(i, j - 1);
    }

    // If right cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      const RealType dx_r         = parameters_.meshsize->getDx(i + 1, j);
      const RealType dx           = parameters_.meshsize->getDx(i, j);
      velocity.getVector(i, j)[1] = -dx / dx_r * velocity.getVector(i + 1, j)[1];
      RealType omega_w = computeW2D(Re, dx_r / 2);
      RealType omega_1 = flowField.getTurbOmega().getScalar(i + 1, j);
      W                = ((dx + dx_r) / dx_r) * (omega_w - omega_1 * (dx / (dx + dx_r)));
      TK               = -dx / dx_r * flowField.getTurbKineticEnergy().getScalar(i + 1, j);
    }
    // Same for left
    if ((obstacle & OBSTACLE_LEFT) == 0) {
      const RealType dx_l         = parameters_.meshsize->getDx(i - 1, j);
      const RealType dx           = parameters_.meshsize->getDx(i, j);
      velocity.getVector(i, j)[1] = -dx / dx_l * velocity.getVector(i - 1, j)[1];
      RealType omega_w = computeW2D(Re, dx_l / 2);
      RealType omega_1 = flowField.getTurbOmega().getScalar(i - 1, j);
      W                = ((dx + dx_l) / dx_l) * (omega_w - omega_1 * (dx / (dx + dx_l)));
      TK               = -dx / dx_l * flowField.getTurbKineticEnergy().getScalar(i - 1, j);
    }

    // Special corner case
    if ((obstacle & OBSTACLE_RIGHT) == 0 && (obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVector(i, j)[0] = 0.0;
      velocity.getVector(i, j)[1] = 0.0;

      RealType& kw = flowField.getTurbKineticEnergy().getScalar(i, j);
      RealType& k1 = flowField.getTurbKineticEnergy().getScalar(i, j + 1);
      RealType& k2 = flowField.getTurbKineticEnergy().getScalar(i + 1, j + 1);
      RealType& k3 = flowField.getTurbKineticEnergy().getScalar(i + 1, j);

      RealType avg = (k1 + k2 + k3) / 3.0;

      kw = -avg;
      k1 = avg;
      k2 = avg;
      k3 = avg;
    }
    // Set normal velocity to zero if right neighbour is not obstacle
    else if ((obstacle & OBSTACLE_RIGHT) == 0) {
      velocity.getVector(i, j)[0] = 0.0;
    }

    // Set normal velocity to zero if top neighbour is not obstacle
    else if ((obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVector(i, j)[1] = 0.0;
    }
  }
}


void Stencils::TurbObstacleStencil::apply([[maybe_unused]] TurbulentFlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k) {}
