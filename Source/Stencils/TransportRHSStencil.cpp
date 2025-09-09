#include "StdAfx.hpp"

#include "TransportRHSStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"


Stencils::TransportRHSStencil::TransportRHSStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::TransportRHSStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const RealType alpha{5.0 / 9.0}, beta{3.0 / 40.0}, beta_star{9.0 / 100.0}, sigma{0.5}, sigma_star{0.5};

  const int    obstacle = flowField.getFlags().getValue(i, j);
  if ((obstacle & OBSTACLE_SELF) == 1) return;


  // Fields to write to
  RealType& RHS_K = flowField.getTurbKineticEnergyRHS().getScalar(i, j);
  RealType& RHS_W = flowField.getTurbOmegaRHS().getScalar(i, j);

  // Needed Fields and values
  VectorField& velocityField  = flowField.getVelocity();
  RealType*    velocityVector = velocityField.getVector(i, j);
  RealType     u              = velocityVector[0];
  RealType     v              = velocityVector[1];

  RealType u_xP = velocityField.getVector(i + 1, j)[0];
  RealType v_xP = velocityField.getVector(i + 1, j)[1];
  RealType u_xN = velocityField.getVector(i - 1, j)[0];
  RealType v_xN = velocityField.getVector(i - 1, j)[1];

  RealType u_yP = velocityField.getVector(i, j + 1)[0];
  RealType v_yP = velocityField.getVector(i, j + 1)[1];
  RealType u_yN = velocityField.getVector(i, j - 1)[0];
  RealType v_yN = velocityField.getVector(i, j - 1)[1];

  RealType u_xP_yP = velocityField.getVector(i + 1, j + 1)[0];
  RealType v_xP_yP = velocityField.getVector(i + 1, j + 1)[1];
  RealType u_xP_yN = velocityField.getVector(i + 1, j - 1)[0];
  RealType v_xP_yN = velocityField.getVector(i + 1, j - 1)[1];

  RealType u_xN_yP = velocityField.getVector(i - 1, j + 1)[0];
  RealType v_xN_yP = velocityField.getVector(i - 1, j + 1)[1];
  RealType u_xN_yN = velocityField.getVector(i - 1, j - 1)[0];
  RealType v_xN_yN = velocityField.getVector(i - 1, j - 1)[1];

  ScalarField& TurbKineticEnergy = flowField.getTurbKineticEnergy();
  RealType     TK                = TurbKineticEnergy.getScalar(i, j);
  RealType     TK_x_next         = TurbKineticEnergy.getScalar(i + 1, j);
  RealType     TK_x_prev         = TurbKineticEnergy.getScalar(i - 1, j);
  RealType     TK_y_next         = TurbKineticEnergy.getScalar(i, j + 1);
  RealType     TK_y_prev         = TurbKineticEnergy.getScalar(i, j - 1);

  ScalarField& TurbOmega    = flowField.getTurbOmega();
  RealType     omega        = TurbOmega.getScalar(i, j);
  RealType     omega_x_next = TurbOmega.getScalar(i + 1, j);
  RealType     omega_x_prev = TurbOmega.getScalar(i - 1, j);
  RealType     omega_y_next = TurbOmega.getScalar(i, j + 1);
  RealType     omega_y_prev = TurbOmega.getScalar(i, j - 1);

  ScalarField& viscosityField = flowField.getTurbViscosity();
  RealType     nu_T           = viscosityField.getScalar(i, j);
  RealType     nu_T_x_next    = viscosityField.getScalar(i + 1, j);
  RealType     nu_T_x_prev    = viscosityField.getScalar(i - 1, j);
  RealType     nu_T_y_next    = viscosityField.getScalar(i, j + 1);
  RealType     nu_T_y_prev    = viscosityField.getScalar(i, j - 1);
  RealType     nu_T_xHP       = 0.5 * (nu_T_x_next + nu_T);
  RealType     nu_T_xHN       = 0.5 * (nu_T_x_prev + nu_T);
  RealType     nu_T_yHP       = 0.5 * (nu_T_y_next + nu_T);
  RealType     nu_T_yHN       = 0.5 * (nu_T_y_prev + nu_T);

  // Meshsizes and Re
  RealType Dx        = parameters_.meshsize->getDx(i, j);
  RealType Dx_x_next = parameters_.meshsize->getDx(i + 1, j);
  RealType Dx_x_prev = parameters_.meshsize->getDx(i - 1, j);
  RealType Dy        = parameters_.meshsize->getDy(i, j);
  RealType Dy_y_next = parameters_.meshsize->getDy(i, j + 1);
  RealType Dy_y_prev = parameters_.meshsize->getDy(i, j - 1);
  RealType Dx_x_h    = 0.5 * (Dx_x_prev + Dx_x_next) + Dx;
  RealType Dy_y_h    = 0.5 * (Dy_y_prev + Dy_y_next) + Dy;
  RealType Re        = parameters_.flow.Re;

  // Gradients and Sij
  RealType dudx = (u - u_xN) / Dx;
  RealType dvdy = (v - v_yN) / Dy;
  RealType dudy = ((u_yP + u_xN_yP) - (u_yN + u_xN_yN)) / (2 * Dy_y_h);
  RealType dvdx = ((v_xP + v_xP_yN) - (v_xN + v_xN_yN)) / (2 * Dx_x_h);

  RealType SQuad = dudx * dudx + 0.5 * (dudy + dvdx) * (dudy + dvdx) + dvdy * dvdy;


  RealType TKRHS_1 = -u * (TK_x_next - TK_x_prev) / Dx_x_h - v * (TK_y_next - TK_y_prev) / Dy_y_h;
  RealType TKRHS_2 = 2 * nu_T * SQuad;
  RealType TKRHS_3 = -beta_star * TK * omega;
  RealType TKRHS_4 = ((1.0 / Re + sigma_star * nu_T_xHP) * (TK_x_next - TK) / (0.5 * (Dx_x_next + Dx))
                      - (1.0 / Re + sigma_star * nu_T_xHN) * (TK - TK_x_prev) / (0.5 * (Dx_x_prev + Dx)))
                     / Dx;
  RealType TKRHS_5 = ((1.0 / Re + sigma_star * nu_T_yHP) * (TK_y_next - TK) / (0.5 * (Dy_y_next + Dy))
                      - (1.0 / Re + sigma_star * nu_T_yHN) * (TK - TK_y_prev) / (0.5 * (Dy_y_prev + Dy)))
                     / Dy;

  RHS_K = TKRHS_1 + TKRHS_2 + TKRHS_3 + TKRHS_4 + TKRHS_5;

  RealType WRHS_1 = -u * (omega_x_next - omega_x_prev) / Dx_x_h - v * (omega_y_next - omega_y_prev) / Dy_y_h;
  RealType WRHS_2 = 2 * alpha * SQuad;
  RealType WRHS_3 = -beta * omega * omega;
  RealType WRHS_4 = ((1.0 / Re + sigma * nu_T_xHP) * (omega_x_next - omega) / (0.5 * (Dx_x_next + Dx))
                     - (1.0 / Re + sigma * nu_T_xHN) * (omega - omega_x_prev) / (0.5 * (Dx_x_prev + Dx)))
                    / Dx;
  RealType WRHS_5 = ((1.0 / Re + sigma * nu_T_yHP) * (omega_y_next - omega) / (0.5 * (Dy_y_next + Dy))
                     - (1.0 / Re + sigma * nu_T_yHN) * (omega - omega_y_prev) / (0.5 * (Dy_y_prev + Dy)))
                    / Dy;

  RHS_W = WRHS_1 + WRHS_2 + WRHS_3 + WRHS_4 + WRHS_5;
}


void Stencils::TransportRHSStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {

  const RealType alpha{5.0 / 9.0}, beta{3.0 / 40.0}, beta_star{9.0 / 100.0}, sigma{0.5}, sigma_star{0.5};

  // Fields to write to
  RealType& RHS_K = flowField.getTurbKineticEnergyRHS().getScalar(i, j, k);
  RealType& RHS_W = flowField.getTurbOmegaRHS().getScalar(i, j, k);

  // Needed Fields and values
  VectorField& velocityField  = flowField.getVelocity();
  RealType*    velocityVector = velocityField.getVector(i, j, k);
  RealType     u              = velocityVector[0];
  RealType     v              = velocityVector[1];
  RealType     w              = velocityVector[2];


  RealType u_xP = velocityField.getVector(i + 1, j, k)[0];
  RealType v_xP = velocityField.getVector(i + 1, j, k)[1];
  RealType w_xP = velocityField.getVector(i + 1, j, k)[2];

  RealType u_xN = velocityField.getVector(i - 1, j, k)[0];
  RealType v_xN = velocityField.getVector(i - 1, j, k)[1];
  RealType w_xN = velocityField.getVector(i - 1, j, k)[2];

  RealType u_yP = velocityField.getVector(i, j + 1, k)[0];
  RealType v_yP = velocityField.getVector(i, j + 1, k)[1];
  RealType w_yP = velocityField.getVector(i, j + 1, k)[2];

  RealType u_yN = velocityField.getVector(i, j - 1, k)[0];
  RealType v_yN = velocityField.getVector(i, j - 1, k)[1];
  RealType w_yN = velocityField.getVector(i, j - 1, k)[2];


  RealType u_zP = velocityField.getVector(i, j, k + 1)[0];
  RealType v_zP = velocityField.getVector(i, j, k + 1)[1];

  RealType w_zN = velocityField.getVector(i, j, k - 1)[2];
  RealType v_zN = velocityField.getVector(i, j, k - 1)[1];
  RealType u_zN = velocityField.getVector(i, j, k - 1)[0];


  RealType u_xP_yP = velocityField.getVector(i + 1, j + 1, k)[0];
  RealType v_xP_yP = velocityField.getVector(i + 1, j + 1, k)[1];

  RealType u_xP_yN = velocityField.getVector(i + 1, j - 1, k)[0];
  RealType v_xP_yN = velocityField.getVector(i + 1, j - 1, k)[1];

  RealType u_xN_zN = velocityField.getVector(i - 1, j, k - 1)[0];
  RealType w_xN_zN = velocityField.getVector(i - 1, j, k - 1)[2];


  RealType u_xN_yP = velocityField.getVector(i - 1, j + 1, k)[0];
  RealType v_xN_yP = velocityField.getVector(i - 1, j + 1, k)[1];

  RealType u_xN_yN = velocityField.getVector(i - 1, j - 1, k)[0];
  RealType v_xN_yN = velocityField.getVector(i - 1, j - 1, k)[1];

  RealType v_yN_zN = velocityField.getVector(i, j - 1, k - 1)[1];
  RealType w_yN_zN = velocityField.getVector(i, j - 1, k - 1)[2];


  RealType v_yN_zP = velocityField.getVector(i, j - 1, k + 1)[1];
  RealType w_yP_zN = velocityField.getVector(i, j + 1, k - 1)[2];
  RealType w_xP_zN = velocityField.getVector(i + 1, j, k - 1)[2];
  RealType u_xN_zP = velocityField.getVector(i - 1, j, k + 1)[0];

  ScalarField& TurbKineticEnergy = flowField.getTurbKineticEnergy();
  RealType     TK                = TurbKineticEnergy.getScalar(i, j, k);
  RealType     TK_x_next         = TurbKineticEnergy.getScalar(i + 1, j, k);
  RealType     TK_x_prev         = TurbKineticEnergy.getScalar(i - 1, j, k);
  RealType     TK_y_next         = TurbKineticEnergy.getScalar(i, j + 1, k);
  RealType     TK_y_prev         = TurbKineticEnergy.getScalar(i, j - 1, k);
  RealType     TK_z_next         = TurbKineticEnergy.getScalar(i, j, k + 1);
  RealType     TK_z_prev         = TurbKineticEnergy.getScalar(i, j, k - 1);

  ScalarField& TurbOmega    = flowField.getTurbOmega();
  RealType     omega        = TurbOmega.getScalar(i, j, k);
  RealType     omega_x_next = TurbOmega.getScalar(i + 1, j, k);
  RealType     omega_x_prev = TurbOmega.getScalar(i - 1, j, k);
  RealType     omega_y_next = TurbOmega.getScalar(i, j + 1, k);
  RealType     omega_y_prev = TurbOmega.getScalar(i, j - 1, k);
  RealType     omega_z_next = TurbOmega.getScalar(i, j, k + 1);
  RealType     omega_z_prev = TurbOmega.getScalar(i, j, k - 1);

  ScalarField& viscosityField = flowField.getTurbViscosity();
  RealType     nu_T           = viscosityField.getScalar(i, j, k);
  RealType     nu_T_x_next    = viscosityField.getScalar(i + 1, j, k);
  RealType     nu_T_x_prev    = viscosityField.getScalar(i - 1, j, k);
  RealType     nu_T_y_next    = viscosityField.getScalar(i, j + 1, k);
  RealType     nu_T_y_prev    = viscosityField.getScalar(i, j - 1, k);
  RealType     nu_T_z_next    = viscosityField.getScalar(i, j, k + 1);
  RealType     nu_T_z_prev    = viscosityField.getScalar(i, j, k - 1);
  RealType     nu_T_xHP       = 0.5 * (nu_T_x_next + nu_T);
  RealType     nu_T_xHN       = 0.5 * (nu_T_x_prev + nu_T);
  RealType     nu_T_yHP       = 0.5 * (nu_T_y_next + nu_T);
  RealType     nu_T_yHN       = 0.5 * (nu_T_y_prev + nu_T);
  RealType     nu_T_zHP       = 0.5 * (nu_T_z_next + nu_T);
  RealType     nu_T_zHN       = 0.5 * (nu_T_z_prev + nu_T);

  // Meshsizes and Re
  RealType Dx        = parameters_.meshsize->getDx(i, j, k);
  RealType Dx_x_next = parameters_.meshsize->getDx(i + 1, j, k);
  RealType Dx_x_prev = parameters_.meshsize->getDx(i - 1, j, k);
  RealType Dy        = parameters_.meshsize->getDy(i, j, k);
  RealType Dy_y_next = parameters_.meshsize->getDy(i, j + 1, k);
  RealType Dy_y_prev = parameters_.meshsize->getDy(i, j - 1, k);
  RealType Dz        = parameters_.meshsize->getDz(i, j, k);
  RealType Dz_z_next = parameters_.meshsize->getDz(i, j, k + 1);
  RealType Dz_z_prev = parameters_.meshsize->getDz(i, j, k - 1);
  RealType Dx_x_h    = 0.5 * (Dx_x_prev + Dx_x_next) + Dx;
  RealType Dy_y_h    = 0.5 * (Dy_y_prev + Dy_y_next) + Dy;
  RealType Dz_z_h    = 0.5 * (Dz_z_prev + Dz_z_next) + Dz;
  RealType Re        = parameters_.flow.Re;


  // Gradients and Sij
  RealType dudx = (u - u_xN) / Dx;
  RealType dvdy = (v - v_yN) / Dy;
  RealType dwdz = (w - w_zN) / Dz;
  RealType dudy = ((u_yP + u_xN_yP) - (u_yN + u_xN_yN)) / (2 * Dy_y_h);
  RealType dudz = ((u_zP + u_xN_zP) - (u_zN + u_xN_zN)) / (2 * Dz_z_h);
  RealType dvdx = ((v_xP + v_xP_yN) - (v_xN + v_xN_yN)) / (2 * Dx_x_h);
  RealType dvdz = ((v_zP + v_yN_zP) - (v_zN + v_yN_zN)) / (2 * Dz_z_h);
  RealType dwdx = ((w_xP + w_xP_zN) - (w_xN + w_xN_zN)) / (2 * Dx_x_h);
  RealType dwdy = ((w_yP + w_yP_zN) - (w_yN + w_yN_zN)) / (2 * Dy_y_h);


  RealType SQuad = dudx * dudx + 0.5 * (dudy + dvdx) * (dudy + dvdx) + dvdy * dvdy + 0.5 * (dudz + dwdx) * (dudz + dwdx) + 0.5 * (dvdz + dwdy) * (dvdz + dwdy) + dwdz * dwdz;


  RealType TKRHS_1 = -u * (TK_x_next - TK_x_prev) / Dx_x_h - v * (TK_y_next - TK_y_prev) / Dy_y_h - w * (TK_z_next - TK_z_prev) / Dz_z_h;
  RealType TKRHS_2 = 2 * nu_T * SQuad;
  RealType TKRHS_3 = -beta_star * TK * omega;
  RealType TKRHS_4 = ((1.0 / Re + sigma_star * nu_T_xHP) * (TK_x_next - TK) / (0.5 * (Dx_x_next + Dx))
                      - (1.0 / Re + sigma_star * nu_T_xHN) * (TK - TK_x_prev) / (0.5 * (Dx_x_prev + Dx)))
                     / Dx;
  RealType TKRHS_5 = ((1.0 / Re + sigma_star * nu_T_yHP) * (TK_y_next - TK) / (0.5 * (Dy_y_next + Dy))
                      - (1.0 / Re + sigma_star * nu_T_yHN) * (TK - TK_y_prev) / (0.5 * (Dy_y_prev + Dy)))
                     / Dy;
  RealType TKRHS_6 = ((1.0 / Re + sigma_star * nu_T_zHP) * (TK_z_next - TK) / (0.5 * (Dz_z_next + Dz))
                      - (1.0 / Re + sigma_star * nu_T_zHN) * (TK - TK_z_prev) / (0.5 * (Dz_z_prev + Dz)))
                     / Dz;

  RHS_K = TKRHS_1 + TKRHS_2 + TKRHS_3 + TKRHS_4 + TKRHS_5 + TKRHS_6;

  RealType WRHS_1 = -u * (omega_x_next - omega_x_prev) / Dx_x_h - v * (omega_y_next - omega_y_prev) / Dy_y_h - w * (omega_z_next - omega_z_prev) / Dz_z_h;
  RealType WRHS_2 = 2 * alpha * SQuad;
  RealType WRHS_3 = -beta * omega * omega;
  RealType WRHS_4 = ((1.0 / Re + sigma * nu_T_xHP) * (omega_x_next - omega) / (0.5 * (Dx_x_next + Dx))
                     - (1.0 / Re + sigma * nu_T_xHN) * (omega - omega_x_prev) / (0.5 * (Dx_x_prev + Dx)))
                    / Dx;
  RealType WRHS_5 = ((1.0 / Re + sigma * nu_T_yHP) * (omega_y_next - omega) / (0.5 * (Dy_y_next + Dy))
                     - (1.0 / Re + sigma * nu_T_yHN) * (omega - omega_y_prev) / (0.5 * (Dy_y_prev + Dy)))
                    / Dy;
  RealType WRHS_6 = ((1.0 / Re + sigma * nu_T_zHP) * (omega_z_next - omega) / (0.5 * (Dz_z_next + Dz))
                     - (1.0 / Re + sigma * nu_T_zHN) * (omega - omega_z_prev) / (0.5 * (Dz_z_prev + Dz)))
                    / Dz;

  RHS_W = WRHS_1 + WRHS_2 + WRHS_3 + WRHS_4 + WRHS_5 + WRHS_6;
}
