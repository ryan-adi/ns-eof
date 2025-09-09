#pragma once

#include "StdAfx.hpp"

#include "Definitions.hpp"
#include "Parameters.hpp"


namespace Stencils {

  // Load the local velocity cube with relevant velocities of the 2D plane
  inline void loadLocalVelocity2D(FlowField& flowField, RealType* const localVelocity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        const RealType* const point                  = flowField.getVelocity().getVector(i + column, j + row);
        localVelocity[39 + 9 * row + 3 * column]     = point[0]; // x-component
        localVelocity[39 + 9 * row + 3 * column + 1] = point[1]; // y-component
      }
    }
  }


  // Load the local velocity cube with surrounding velocities
  inline void loadLocalVelocity3D(FlowField& flowField, RealType* const localVelocity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          const RealType* const point                               = flowField.getVelocity().getVector(i + column, j + row, k + layer);
          localVelocity[39 + 27 * layer + 9 * row + 3 * column]     = point[0]; // x-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 1] = point[1]; // y-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 2] = point[2]; // z-component
        }
      }
    }
  }

  // Load local meshsize for 2D -> same as loadLocalVelocity2D, but invoking call to meshsize-ptr
  inline void loadLocalMeshsize2D(const Parameters& parameters, RealType* const localMeshsize, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localMeshsize[39 + 9 * row + 3 * column]     = parameters.meshsize->getDx(i + column, j + row);
        localMeshsize[39 + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(i + column, j + row);
      }
    }
  }

  // Load local meshsize for 3D
  inline void loadLocalMeshsize3D(const Parameters& parameters, RealType* const localMeshsize, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column]     = parameters.meshsize->getDx(i + column, j + row, k + layer);
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(i + column, j + row, k + layer);
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 2] = parameters.meshsize->getDz(i + column, j + row, k + layer);
        }
      }
    }
  }

  inline void loadLocalViscosity2D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localViscosity[39 + 9 * row + 3 * column] = flowField.getTurbViscosity().getScalar(i + column, j + row);
      }
    }
  }
  // inline void loadLocalTurbKineticEnergy2D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j) {
  //   for (int row = -1; row <= 1; row++) {
  //     for (int column = -1; column <= 1; column++) {
  //       LocalTurbKineticEnergy[39 + 9 * row + 3 * column] = flowField.getTurbKineticEnergy().getScalar(i + column, j + row);
  //     }
  //   }
  // }

  // inline void loadLocalTurbKineticEnergy3D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j, int k) {
  //   for (int layer = -1; layer <= 1; layer++) {
  //     for (int row = -1; row <= 1; row++) {
  //       for (int column = -1; column <= 1; column++) {
  //         LocalTurbKineticEnergy[39 + 27 * layer + 9 * row + 3 * column] = flowField.getTurbKineticEnergy().getScalar(i + column, j + row, k + layer);
  //       }
  //     }
  //   }
  // }


  inline void loadLocalViscosity3D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localViscosity[39 + 27 * layer + 9 * row + 3 * column] = flowField.getTurbViscosity().getScalar(i + column, j + row, k + layer);
        }
      }
    }
  }

  // Maps an index and a component to the corresponding value in the cube.
  inline int mapd(int i, int j, int k, int component) { return 39 + 27 * k + 9 * j + 3 * i + component; }

  inline RealType F1(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate F1 in the cell center by a central difference.
    const int x0_0  = mapd(0, 0, 0, 0);
    const int xP1_0 = mapd(1, 0, 0, 0);
    const int xM1_0 = mapd(-1, 0, 0, 0);

    // no Re part gets deleted by differentiation
    return ((tv[xP1_0] + 1.0 / Re) * (lv[xP1_0] - lv[x0_0]) / lm[xP1_0] - (tv[x0_0] + 1.0 / Re) * (lv[x0_0] - lv[xM1_0]) / lm[x0_0]) / (0.5 * (lm[x0_0] + lm[xP1_0]));
  }


  inline RealType G2(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate G2 in the cell center by a central difference
    const int y0_1  = mapd(0, 0, 0, 1);
    const int yP1_1 = mapd(0, 1, 0, 1);
    const int yM1_1 = mapd(0, -1, 0, 1);

    const int y0_0  = mapd(0, 0, 0, 0); // viscostiy is stored in 0
    const int yP1_0 = mapd(0, 1, 0, 0); // viscostiy is stored in 0


    // no Re part gets deleted by differentiation
    return ((tv[yP1_0] + 1.0 / Re) * (lv[yP1_1] - lv[y0_1]) / lm[yP1_1] - (tv[y0_0] + 1.0 / Re) * (lv[y0_1] - lv[yM1_1]) / lm[y0_1]) / (0.5 * (lm[yP1_1] + lm[y0_1]));
  }

  inline RealType H3(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate H3 in the cell center by a central difference
    const int z0_2  = mapd(0, 0, 0, 2);
    const int zP1_2 = mapd(0, 0, 1, 2);
    const int zM1_2 = mapd(0, 0, -1, 2);

    const int z0_0  = mapd(0, 0, 0, 0); // viscostiy is stored in 0s
    const int zP1_0 = mapd(0, 0, 1, 0); // viscostiy is stored in 0

    // no Re part gets deleted by differentiation
    return ((tv[zP1_0] + 1.0 / Re) * (lv[zP1_2] - lv[z0_2]) / lm[zP1_2] - (tv[z0_0] + 1.0 / Re) * (lv[z0_2] - lv[zM1_2]) / lm[z0_2]) / (0.5 * (lm[zP1_2] + lm[z0_2]));
  }


  inline RealType F2(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate F2 in the cell center by a central difference
    const int x0_y0_0   = mapd(0, 0, 0, 0);
    const int xP1_y0_0  = mapd(1, 0, 0, 0);
    const int x0_yM1_0  = mapd(0, -1, 0, 0);
    const int xP1_yM1_0 = mapd(1, -1, 0, 0);
    const int x0_yP1_0  = mapd(0, 1, 0, 0);
    const int xP1_yP1_0 = mapd(1, 1, 0, 0);

    const int x0_yP1_1  = mapd(0, 1, 0, 1);
    const int x0_y0_1   = mapd(0, 0, 0, 1);
    const int xP1_y0_1  = mapd(1, 0, 0, 1);
    const int x0_yM1_1  = mapd(0, -1, 0, 1);
    const int xP1_yM1_1 = mapd(1, -1, 0, 1);

    const RealType atv_minus = (tv[x0_y0_0] + tv[xP1_y0_0] + tv[x0_yM1_0] + tv[xP1_yM1_0]) / 4.0 + 1.0 / Re;
    const RealType atv_plus  = (tv[x0_y0_0] + tv[xP1_y0_0] + tv[x0_yP1_0] + tv[xP1_yP1_0]) / 4.0 + 1.0 / Re;

    const RealType lm_avg1 = 0.5 * (lm[x0_yP1_1] + lm[x0_y0_1]);
    const RealType lm_avg2 = 0.5 * (lm[xP1_y0_0] + lm[x0_y0_0]);
    const RealType lm_avg3 = 0.5 * (lm[x0_y0_1] + lm[x0_yM1_1]);
    const RealType lm_avg4 = 0.5 * (lm[xP1_yM1_0] + lm[x0_yM1_0]);


    return ((atv_plus) * ((lv[x0_yP1_0] - lv[x0_y0_0]) / lm_avg1 + (lv[xP1_y0_1] - lv[x0_y0_1]) / lm_avg2)
            - (atv_minus) * ((lv[x0_y0_0] - lv[x0_yM1_0]) / lm_avg3 + (lv[xP1_yM1_1] - lv[x0_yM1_1]) / lm_avg4))
           / (lm[x0_y0_1]);
  }

  inline RealType G1(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate G1 in the cell center by a central difference
    const int x0_y0_0   = mapd(0, 0, 0, 0);
    const int xP1_y0_0  = mapd(1, 0, 0, 0);
    const int xM1_y0_0  = mapd(-1, 0, 0, 0);
    const int xM1_yP1_0 = mapd(-1, 1, 0, 0);
    const int x0_yP1_0  = mapd(0, 1, 0, 0);
    const int xP1_yP1_0 = mapd(1, 1, 0, 0);

    const int x0_y0_1   = mapd(0, 0, 0, 1);
    const int x0_yP1_1  = mapd(0, 1, 0, 1);
    const int xP1_y0_1  = mapd(1, 0, 0, 1);
    const int xM1_y0_1  = mapd(-1, 0, 0, 1);
    const int xM1_yP1_1 = mapd(-1, 1, 0, 1);

    const RealType atv_minus = (tv[x0_y0_0] + tv[x0_yP1_0] + tv[xM1_y0_0] + tv[xM1_yP1_0]) / 4.0 + 1.0 / Re;
    const RealType atv_plus  = (tv[x0_y0_0] + tv[xP1_y0_0] + tv[x0_yP1_0] + tv[xP1_yP1_0]) / 4.0 + 1.0 / Re;

    const RealType lm_avg1 = 0.5 * (lm[x0_yP1_1] + lm[x0_y0_1]);
    const RealType lm_avg2 = 0.5 * (lm[xP1_y0_0] + lm[x0_y0_0]);
    const RealType lm_avg3 = 0.5 * (lm[xM1_yP1_1] + lm[xM1_y0_1]);
    const RealType lm_avg4 = 0.5 * (lm[x0_y0_0] + lm[xM1_y0_0]);

    return ((atv_plus) * ((lv[x0_yP1_0] - lv[x0_y0_0]) / lm_avg1 + (lv[xP1_y0_1] - lv[x0_y0_1]) / lm_avg2)
            - (atv_minus) * ((lv[xM1_yP1_0] - lv[xM1_y0_0]) / lm_avg3 + (lv[x0_y0_1] - lv[xM1_y0_1]) / lm_avg4))
           / (lm[x0_y0_0]);
  }


  inline RealType H1(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate dudx in the cell center by a central difference
    const int x0_z0_0     = mapd(0, 0, 0, 0);
    const int x_P1_z0_0   = mapd(1, 0, 0, 0);
    const int x_M1_z0_0   = mapd(-1, 0, 0, 0);
    const int x_M1_z_P1_0 = mapd(-1, 0, 1, 0);
    const int x0_z_P1_0   = mapd(0, 0, 1, 0);
    const int x_P1_z_P1_0 = mapd(1, 0, 1, 0);


    const int x0_z0_2     = mapd(0, 0, 0, 2);
    const int x_P1_z0_2   = mapd(1, 0, 0, 2);
    const int x_M1_z0_2   = mapd(-1, 0, 0, 2);
    const int x0_z_P1_2   = mapd(0, 0, 1, 2);
    const int x_M1_z_P1_2 = mapd(-1, 0, 1, 2);

    const RealType atv_minus = (tv[x0_z0_0] + tv[x_M1_z0_0] + tv[x_M1_z_P1_0] + tv[x0_z_P1_0]) / 4.0 + 1.0 / Re;
    const RealType atv_plus  = (tv[x0_z0_0] + tv[x_P1_z0_0] + tv[x0_z_P1_0] + tv[x_P1_z_P1_0]) / 4.0 + 1.0 / Re;

    const RealType lm_avg1 = 0.5 * (lm[x_P1_z0_0] + lm[x0_z0_0]);
    const RealType lm_avg2 = 0.5 * (lm[x0_z_P1_2] + lm[x0_z0_2]);
    const RealType lm_avg3 = 0.5 * (lm[x0_z0_0] + lm[x_M1_z0_0]);
    const RealType lm_avg4 = 0.5 * (lm[x_M1_z_P1_2] + lm[x_M1_z0_2]);

    return ((atv_plus) * ((lv[x_P1_z0_2] - lv[x0_z0_2]) / lm_avg1 + (lv[x0_z_P1_0] - lv[x0_z0_0]) / lm_avg2)
            - (atv_minus) * ((lv[x0_z0_2] - lv[x_M1_z0_2]) / lm_avg3 + (lv[x_M1_z_P1_0] - lv[x_M1_z0_0]) / lm_avg4))
           / (lm[x0_z0_0]);
  }


  inline RealType F3(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate dudx in the cell center by a central difference
    const int x0_z0_0     = mapd(0, 0, 0, 0);
    const int x_P1_z0_0   = mapd(1, 0, 0, 0);
    const int x0_z_M1_0   = mapd(0, 0, -1, 0);
    const int x_P1_z_M1_0 = mapd(1, 0, -1, 0);
    const int x0_z_P1_0   = mapd(0, 0, 1, 0);
    const int x_P1_z_P1_0 = mapd(1, 0, 1, 0);


    const int x0_z0_2   = mapd(0, 0, 0, 2);
    const int x0_z_P1_2 = mapd(0, 0, 1, 2);

    const int x_P1_z0_2   = mapd(1, 0, 0, 2);
    const int x0_z_M1_2   = mapd(0, 0, -1, 2);
    const int x_P1_z_M1_2 = mapd(1, 0, -1, 2);


    const RealType atv_minus = (tv[x0_z0_0] + tv[x_P1_z0_0] + tv[x_P1_z_M1_0] + tv[x0_z_M1_0]) / 4.0 + 1.0 / Re;
    const RealType atv_plus  = (tv[x0_z0_0] + tv[x_P1_z0_0] + tv[x0_z_P1_0] + tv[x_P1_z_P1_0]) / 4.0 + 1.0 / Re;

    const RealType lm_avg1 = 0.5 * (lm[x0_z_P1_2] + lm[x0_z0_2]);
    const RealType lm_avg2 = 0.5 * (lm[x_P1_z0_0] + lm[x0_z0_0]);
    const RealType lm_avg3 = 0.5 * (lm[x0_z0_2] + lm[x0_z_M1_2]);
    const RealType lm_avg4 = 0.5 * (lm[x_P1_z_M1_0] + lm[x0_z_M1_0]);

    return ((atv_plus) * ((lv[x0_z_P1_0] - lv[x0_z0_0]) / lm_avg1 + (lv[x_P1_z0_2] - lv[x0_z0_2]) / lm_avg2)
            - (atv_minus) * ((lv[x0_z0_0] - lv[x0_z_M1_0]) / lm_avg3 + (lv[x_P1_z_M1_2] - lv[x0_z_M1_2]) / lm_avg4))
           / (lm[x0_z0_2]);
  }


  inline RealType G3(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate dudx in the cell center by a central difference
    const int y0_z0_1     = mapd(0, 0, 0, 1);
    const int y_P1_z0_1   = mapd(0, 1, 0, 1);
    const int y0_z_M1_1   = mapd(0, 0, -1, 1);
    const int y_P1_z_M1_1 = mapd(0, 1, -1, 1);
    const int y0_z_P1_1   = mapd(0, 0, 1, 1);

    const int y0_z0_2     = mapd(0, 0, 0, 2);
    const int y_P1_z0_2   = mapd(0, 1, 0, 2);
    const int y0_z_M1_2   = mapd(0, 0, -1, 2);
    const int y0_z_P1_2   = mapd(0, 0, 1, 2);
    const int y_P1_z_M1_2 = mapd(0, 1, -1, 2);

    const int index0 = mapd(0, 0, 0, 0);
    const int index1 = mapd(0, 1, 0, 0);
    const int index2 = mapd(0, 0, -1, 0);
    const int index3 = mapd(0, 1, -1, 0);
    const int index4 = mapd(0, 0, 1, 0);
    const int index5 = mapd(0, 1, 1, 0);

    const RealType lm_avg1 = 0.5 * (lm[y0_z_P1_2] + lm[y0_z0_2]);
    const RealType lm_avg2 = 0.5 * (lm[y_P1_z0_1] + lm[y0_z0_1]);
    const RealType lm_avg3 = 0.5 * (lm[y0_z0_2] + lm[y0_z_M1_2]);
    const RealType lm_avg4 = 0.5 * (lm[y_P1_z_M1_1] + lm[y0_z_M1_1]);

    const RealType atv_minus = (tv[index0] + tv[index1] + tv[index2] + tv[index3]) / 4.0 + 1.0 / Re;
    const RealType atv_plus  = (tv[index0] + tv[index1] + tv[index4] + tv[index5]) / 4.0 + 1.0 / Re;

    return ((atv_plus) * ((lv[y0_z_P1_1] - lv[y0_z0_1]) / lm_avg1 + (lv[y_P1_z0_2] - lv[y0_z0_2]) / lm_avg2)
            - (atv_minus) * ((lv[y0_z0_1] - lv[y0_z_P1_1]) / lm_avg3 + (lv[y_P1_z_M1_2] - lv[y0_z_M1_2]) / lm_avg4))
           / (lm[y0_z0_2]);
  }


  inline RealType H2(const RealType* const tv, const RealType* const lm, const RealType* const lv, const RealType Re) {
    // Evaluate dudx in the cell center by a central difference
    const int y0_z0_2     = mapd(0, 0, 0, 2);
    const int y0_z_P1_2   = mapd(0, 0, 1, 2);
    const int y_M1_z_P1_2 = mapd(0, -1, 1, 2);
    // const int y_P1_z_M1_2 = mapd(0, -1, 1, 2);
    const int y_P1_z0_2 = mapd(0, 1, 0, 2);
    const int y_M1_z0_2 = mapd(0, -1, 0, 2);

    const int y0_z0_1     = mapd(0, 0, 0, 1);
    const int y_P1_z0_1   = mapd(0, 1, 0, 1);
    const int y_M1_z0_1   = mapd(0, -1, 0, 1);
    const int y_M1_z_P1_1 = mapd(0, -1, 1, 1);
    const int y0_z_P1_1   = mapd(0, 0, 1, 1);


    const int index0 = mapd(0, 0, 0, 0);
    const int index1 = mapd(0, 1, 0, 0);
    const int index2 = mapd(0, 0, -1, 0);
    const int index3 = mapd(0, 1, -1, 0);
    const int index4 = mapd(0, 0, 1, 0);
    const int index5 = mapd(0, 1, 1, 0);

    const RealType lm_avg1 = 0.5 * (lm[y_P1_z0_1] + lm[y0_z0_1]);
    const RealType lm_avg2 = 0.5 * (lm[y0_z_P1_2] + lm[y0_z0_2]);
    const RealType lm_avg3 = 0.5 * (lm[y0_z0_1] + lm[y_M1_z_P1_1]);
    const RealType lm_avg4 = 0.5 * (lm[y_M1_z_P1_2] + lm[y_M1_z0_2]);

    const RealType atv_minus = (tv[index0] + tv[index1] + tv[index2] + tv[index3]) / 4.0 + 1.0 / Re;
    const RealType atv_plus  = (tv[index0] + tv[index1] + tv[index4] + tv[index5]) / 4.0 + 1.0 / Re;

    return ((atv_plus) * ((lv[y_P1_z0_2] - lv[y0_z0_2]) / lm_avg1 + (lv[y0_z_P1_1] - lv[y0_z0_1]) / lm_avg2)
            - (atv_minus) * ((lv[y0_z0_2] - lv[y_M1_z_P1_2]) / lm_avg3 + (lv[y_M1_z_P1_1] - lv[y_M1_z0_1]) / lm_avg4))
           / (lm[y0_z0_1]);
  }


  // Derivative functions. They are applied to a cube of 3x3x3 cells. lv stands for the local velocity, lm represents
  // the local mesh sizes dudx <-> first derivative of u-component of velocity field w.r.t. x-direction.
  inline RealType dudx(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int index0 = mapd(1, 0, 0, 0);
    const int index1 = mapd(-1, 0, 0, 0);
    return (lv[index0] - lv[index1]) / (2*lm[index0]);
  }
  inline RealType dvdy(const RealType* const lv, const RealType* const lm) {
    const int index0 = mapd(0, 1, 0, 1);
    const int index1 = mapd(0, -1, 0, 1);
    return (lv[index0] - lv[index1]) / (2*lm[index0]);
  }

  inline RealType dwdz(const RealType* const lv, const RealType* const lm) {
    const int index0 = mapd(0, 0, 1, 2);
    const int index1 = mapd(0, 0, -1, 2);
    return (lv[index0] - lv[index1]) / (2*lm[index0]);
  }

  // Second derivative of u-component w.r.t. x-direction, evaluated at the location of the u-component.
  inline RealType d2udx2(const RealType* const lv, const RealType* const lm) {
    // Evaluate the second derivative at the location of the u-component of the velocity field;
    // we therefore use the two neighbouring u-components and assume arbitrary mesh sizes in both
    // directions -> the formula arises from a straight-forward taylor expansion
    // -> for equal meshsizes, we obtain the usual [1 -2 1]-like stencil.
    const int indexM1 = mapd(-1, 0, 0, 0);
    const int index0  = mapd(0, 0, 0, 0);
    const int indexP1 = mapd(1, 0, 0, 0);

    const RealType dx0   = lm[index0];
    const RealType dx1   = lm[indexP1];
    const RealType dxSum = dx0 + dx1;
    return 2.0 * (lv[indexP1] / (dx1 * dxSum) - lv[index0] / (dx1 * dx0) + lv[indexM1] / (dx0 * dxSum));
  }

  inline RealType d2udy2(const RealType* const lv, const RealType* const lm) {
    // Average mesh sizes, since the component u is located in the middle of the cell's face.
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);
    const RealType dySum = dy0 + dy1;
    return 2.0 * (lv[mapd(0, 1, 0, 0)] / (dy1 * dySum) - lv[mapd(0, 0, 0, 0)] / (dy1 * dy0) + lv[mapd(0, -1, 0, 0)] / (dy0 * dySum));
  }

  inline RealType d2udz2(const RealType* const lv, const RealType* const lm) {
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);
    const RealType dzSum = dz0 + dz1;
    return 2.0 * (lv[mapd(0, 0, 1, 0)] / (dz1 * dzSum) - lv[mapd(0, 0, 0, 0)] / (dz1 * dz0) + lv[mapd(0, 0, -1, 0)] / (dz0 * dzSum));
  }

  // Second derivative of the v-component, evaluated at the location of the v-component.
  inline RealType d2vdx2(const RealType* const lv, const RealType* const lm) {
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);
    const RealType dxSum = dx0 + dx1;
    return 2.0 * (lv[mapd(1, 0, 0, 1)] / (dx1 * dxSum) - lv[mapd(0, 0, 0, 1)] / (dx1 * dx0) + lv[mapd(-1, 0, 0, 1)] / (dx0 * dxSum));
  }

  inline RealType d2vdy2(const RealType* const lv, const RealType* const lm) {
    const int indexM1 = mapd(0, -1, 0, 1);
    const int index0  = mapd(0, 0, 0, 1);
    const int indexP1 = mapd(0, 1, 0, 1);

    const RealType dy0   = lm[index0];
    const RealType dy1   = lm[indexP1];
    const RealType dySum = dy0 + dy1;
    return 2.0 * (lv[indexP1] / (dy1 * dySum) - lv[index0] / (dy1 * dy0) + lv[indexM1] / (dy0 * dySum));
  }

  inline RealType d2vdz2(const RealType* const lv, const RealType* const lm) {
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);
    const RealType dzSum = dz0 + dz1;
    return 2.0 * (lv[mapd(0, 0, 1, 1)] / (dz1 * dzSum) - lv[mapd(0, 0, 0, 1)] / (dz1 * dz0) + lv[mapd(0, 0, -1, 1)] / (dz0 * dzSum));
  }

  // Second derivative of the w-component, evaluated at the location of the w-component.
  inline RealType d2wdx2(const RealType* const lv, const RealType* const lm) {
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);
    const RealType dxSum = dx0 + dx1;
    return 2.0 * (lv[mapd(1, 0, 0, 2)] / (dx1 * dxSum) - lv[mapd(0, 0, 0, 2)] / (dx1 * dx0) + lv[mapd(-1, 0, 0, 2)] / (dx0 * dxSum));
  }

  inline RealType d2wdy2(const RealType* const lv, const RealType* const lm) {
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);
    const RealType dySum = dy0 + dy1;
    return 2.0 * (lv[mapd(0, 1, 0, 2)] / (dy1 * dySum) - lv[mapd(0, 0, 0, 2)] / (dy1 * dy0) + lv[mapd(0, -1, 0, 2)] / (dy0 * dySum));
  }

  inline RealType d2wdz2(const RealType* const lv, const RealType* const lm) {
    const int index_M1 = mapd(0, 0, -1, 2);
    const int index_0  = mapd(0, 0, 0, 2);
    const int index_P1 = mapd(0, 0, 1, 2);

    const RealType dz0   = lm[index_0];
    const RealType dz1   = lm[index_P1];
    const RealType dzSum = dz0 + dz1;
    return 2.0 * (lv[index_P1] / (dz1 * dzSum) - lv[index_0] / (dz1 * dz0) + lv[index_M1] / (dz0 * dzSum));
  }

  // First derivative of product (u*v), evaluated at the location of the v-component.
  inline RealType duvdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 1)] + lv[mapd(0, 0, 0, 1)])))
        + parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(1, 0, 0, 1)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)])))
        ) / lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of center u-value from upper edge of cell
    const RealType hyLong  = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 1, 0, 0)];
    const RealType vM10 = lv[mapd(-1, 0, 0, 1)];

    // This a central difference expression for the first-derivative. We therefore linearly interpolate u*v onto the
    // surface of the current cell (in 2D: upper left and upper right corner) and then take the central difference.
    const RealType secondOrder = (((hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01) * ((hxLong1 - hxShort) / hxLong1 * v00 + hxShort / hxLong1 * v10)
                                  - ((hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11) * ((hxLong0 - hxShort) / hxLong0 * v00 + hxShort / hxLong0 * vM10))
                                 / (2.0 * hxShort);

    // This is a forward-difference in donor-cell style. We apply donor cell and again interpolate the velocity values
    // (u-comp.) onto the surface of the cell. We then apply the standard donor cell scheme. This will, however, result
    // in non-equal mesh spacing evaluations (in case of stretched meshes).
    const RealType kr = (hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01;
    const RealType kl = (hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11;

    const RealType firstOrder = 1.0 / (4.0 * hxShort) * (kr * (v00 + v10) - kl * (vM10 + v00) + fabs(kr) * (v00 - v10) - fabs(kl) * (vM10 - v00));

    // Return linear combination of central and donor-cell difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for u*v at location of u-component. For details on implementation, see duvdx.
  inline RealType duvdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
            (lv[mapd(0, -1, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 1, 0, 0)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
                (lv[mapd(0, -1, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of center u-value from upper edge of cell
    const RealType hxLong  = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];

    const RealType v0M1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1M1 = lv[mapd(1, -1, 0, 1)];
    const RealType u0M1 = lv[mapd(0, -1, 0, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10) * ((hyLong1 - hyShort) / hyLong1 * u00 + hyShort / hyLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1) * ((hyLong0 - hyShort) / hyLong0 * u00 + hyShort / hyLong0 * u0M1))
                                 / (2.0 * hyShort);

    const RealType kr = (hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10;
    const RealType kl = (hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1;

    const RealType firstOrder = 1.0 / (4.0 * hyShort) * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. x for u*w at location of w-component. For details on implementation, see duvdx.
  inline RealType duwdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
            (lv[mapd(-1, 0, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(1, 0, 0, 2)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
                (lv[mapd(-1, 0, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of center u-value from upper edge of cell
    const RealType hzLong  = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 0, 1, 0)];
    const RealType wM10 = lv[mapd(-1, 0, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01) * ((hxLong1 - hxShort) / hxLong1 * w00 + hxShort / hxLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11) * ((hxLong0 - hxShort) / hxLong0 * w00 + hxShort / hxLong0 * wM10))
                                 / (2.0 * hxShort);

    const RealType kr = (hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01;
    const RealType kl = (hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11;

    const RealType firstOrder = 1.0 / (4.0 * hxShort) * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for u*w at location of u-component. For details on implementation, see duvdx.
  inline RealType duwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
            (lv[mapd(0, 0, -1, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 0, 1, 0)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
                (lv[mapd(0, 0, -1, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];                           // Distance of center u-value from upper edge of cell
    const RealType hxLong  = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(1, 0, -1, 2)];
    const RealType u0M1 = lv[mapd(0, 0, -1, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10) * ((hzLong1 - hzShort) / hzLong1 * u00 + hzShort / hzLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1) * ((hzLong0 - hzShort) / hzLong0 * u00 + hzShort / hzLong0 * u0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10;
    const RealType kl = (hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1;

    const RealType firstOrder = 1.0 / (4.0 * hzShort) * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdz");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for v*w at location of w-component. For details on implementation, see duvdx.
  inline RealType dvwdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
            (lv[mapd(0, -1, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 1, 0, 2)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
                (lv[mapd(0, -1, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of center u-value from upper edge of cell
    const RealType hzLong  = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];

    const RealType vM10 = lv[mapd(0, -1, 0, 1)];
    const RealType vM11 = lv[mapd(0, -1, 1, 1)];
    const RealType wM10 = lv[mapd(0, -1, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01) * ((hyLong1 - hyShort) / hyLong1 * w00 + hyShort / hyLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11) * ((hyLong0 - hyShort) / hyLong0 * w00 + hyShort / hyLong0 * wM10))
                                 / (2.0 * hyShort);

    const RealType kr = (hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01;
    const RealType kl = (hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11;

    const RealType firstOrder = 1.0 / (4.0 * hyShort) * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in dvwdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for v*w at location of v-component. For details on implementation, see duvdx.
  inline RealType dvwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
            (lv[mapd(0, 0, -1, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 0, 1, 1)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
                (lv[mapd(0, 0, -1, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];                           // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];                           // Distance of center u-value from upper edge of cell
    const RealType hyLong  = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(0, 1, -1, 2)];
    const RealType v0M1 = lv[mapd(0, 0, -1, 1)];

    const RealType secondOrder = (((hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10) * ((hzLong1 - hzShort) / hzLong1 * v00 + hzShort / hzLong1 * v01)
                                  - ((hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1) * ((hzLong0 - hzShort) / hzLong0 * v00 + hzShort / hzLong0 * v0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10;
    const RealType kl = (hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1;

    const RealType firstOrder = 1.0 / (4.0 * hzShort) * (kr * (v00 + v01) - kl * (v0M1 + v00) + fabs(kr) * (v00 - v01) - fabs(kl) * (v0M1 - v00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dvwdz");
    }
#endif

    return tmp2;
  }

  // First derivative of u*u w.r.t. x, evaluated at location of u-component.
  inline RealType du2dx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(1, 0, 0, 0)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType dxShort = 0.5 * lm[mapd(0, 0, 0, 0)];
    // const RealType dxLong0 = 0.5*(lm[mapd(-1,0,0,0)] + lm[mapd(0,0,0,0)]);
    const RealType dxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);

    const RealType u0  = lv[mapd(0, 0, 0, 0)];
    const RealType uM1 = lv[mapd(-1, 0, 0, 0)];
    const RealType u1  = lv[mapd(1, 0, 0, 0)];

    // const RealType kr = (dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1;
    // const RealType kl = (dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1;
    const RealType kr = (u0 + u1) / 2;
    const RealType kl = (u0 + uM1) / 2;

    // Central difference expression which is second-order accurate for uniform meshes. We interpolate u half-way
    // between neighboured u-component values and afterwards build the central difference for u*u.

    /*const RealType secondOrder = (((dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1) * ((dxLong1 - dxShort)
       / dxLong1 * u0 + dxShort / dxLong1 * u1)
        - ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1) * ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort
       / dxLong0 * uM1) ) / (2.0 * dxShort);*/

    const RealType secondOrder = ((u0 + u1) * (u0 + u1) - (u0 + uM1) * (u0 + uM1)) / (4 * dxLong1);

    // Donor-cell like derivative expression. We evaluate u half-way between neighboured u-components and use this as a
    // prediction of the transport direction.
    const RealType firstOrder = 1.0 / (4.0 * dxShort) * (kr * (u0 + u1) - kl * (uM1 + u0) + fabs(kr) * (u0 - u1) - fabs(kl) * (uM1 - u0));

    // Return linear combination of central- and upwind difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in du2dx");
    }
#endif

    return tmp2;
  }

  // First derivative of v*v w.r.t. y, evaluated at location of v-component. For details, see du2dx.
  inline RealType dv2dy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
            (lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 1, 0, 1)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
                (lv[mapd(0, -1, 0, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType dyShort = 0.5 * lm[mapd(0, 0, 0, 1)];
    // const RealType dyLong0 = 0.5*(lm[mapd(0,-1,0,1)] + lm[mapd(0,0,0,1)]);
    const RealType dyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);

    const RealType v0  = lv[mapd(0, 0, 0, 1)];
    const RealType vM1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1  = lv[mapd(0, 1, 0, 1)];

    // const RealType kr = (dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1;
    // const RealType kl = (dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1;
    const RealType kr = (v0 + v1) / 2;
    const RealType kl = (v0 + vM1) / 2;

    /*const RealType secondOrder = (((dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1) * ((dyLong1 - dyShort)
       / dyLong1 * v0 + dyShort / dyLong1 * v1)
        - ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1) * ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort
       / dyLong0 * vM1) ) / (2.0 * dyShort);*/

    const RealType secondOrder = ((v0 + v1) * (v0 + v1) - (v0 + vM1) * (v0 + vM1)) / (4 * dyLong1);

    const RealType firstOrder = 1.0 / (4.0 * dyShort) * (kr * (v0 + v1) - kl * (vM1 + v0) + fabs(kr) * (v0 - v1) - fabs(kl) * (vM1 - v0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dv2dy");
    }
#endif

    return tmp2;
  }

  // First derivative of w*w w.r.t. z, evaluated at location of w-component. For details, see du2dx.
  inline RealType dw2dz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
            (lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 0, 1, 2)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
                (lv[mapd(0, 0, -1, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType dzShort = 0.5 * lm[mapd(0, 0, 0, 2)];
    // const RealType dzLong0 = 0.5 * (lm[mapd(0, 0, -1, 2)] + lm[mapd(0, 0, 0, 2)]);
    const RealType dzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);

    const RealType w0  = lv[mapd(0, 0, 0, 2)];
    const RealType wM1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1  = lv[mapd(0, 0, 1, 2)];

    // const RealType kr = (dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1;
    // const RealType kl = (dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1;
    const RealType kr = (w0 + w1) / 2;
    const RealType kl = (w0 + wM1) / 2;

    /*const RealType secondOrder = (((dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1) * ((dzLong1 - dzShort)
       / dzLong1 * w0 + dzShort / dzLong1 * w1)
        - ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1) * ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort
       / dzLong0 * wM1) ) / (2.0 * dzShort);*/

    const RealType secondOrder = ((w0 + w1) * (w0 + w1) - (w0 + wM1) * (w0 + wM1)) / (4 * dzLong1);

    const RealType firstOrder = 1.0 / (4.0 * dzShort) * (kr * (w0 + w1) - kl * (wM1 + w0) + fabs(kr) * (w0 - w1) - fabs(kl) * (wM1 - w0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dw2dz");
    }
#endif

    return tmp2;
  }


  inline RealType computeF2D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {

    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (1 / parameters.flow.Re * (d2udx2(localVelocity, localMeshsize)
            + d2udy2(localVelocity, localMeshsize)) - du2dx(localVelocity, parameters, localMeshsize)
            - duvdy(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeG2D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (1 / parameters.flow.Re * (d2vdx2(localVelocity, localMeshsize)
            + d2vdy2(localVelocity, localMeshsize)) - duvdx(localVelocity, parameters, localMeshsize)
            - dv2dy(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computeF3D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (1 / parameters.flow.Re * (d2udx2(localVelocity, localMeshsize)
            + d2udy2(localVelocity, localMeshsize) + d2udz2(localVelocity, localMeshsize))
            - du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize)
            - duwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeG3D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (1 / parameters.flow.Re * (d2vdx2(localVelocity, localMeshsize)
            + d2vdy2(localVelocity, localMeshsize) + d2vdz2(localVelocity, localMeshsize))
            - dv2dy(localVelocity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize)
            - dvwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computeH3D(const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt) {
    return localVelocity[mapd(0, 0, 0, 2)]
        + dt * (1 / parameters.flow.Re * (d2wdx2(localVelocity, localMeshsize)
            + d2wdy2(localVelocity, localMeshsize) + d2wdz2(localVelocity, localMeshsize))
            - dw2dz(localVelocity, parameters, localMeshsize) - duwdx(localVelocity, parameters, localMeshsize)
            - dvwdy(localVelocity, parameters, localMeshsize) + parameters.environment.gz);
  }

  ///////////////////////////////////////////////////////RYAN////////////////////////////////////////////////////////
  // WS2 new local derivatives for calculation of S_ij * S_ij in turbulent viscosity equation (s. Griebel)
  inline RealType dudy(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudy in the cell center by a central difference
    const int index0 = mapd(0, 1, 0, 0);
    const int index1 = mapd(-1, 1, 0, 0);
    const int index2 = mapd(0, -1, 0, 0);
    const int index3 = mapd(-1, -1, 0, 0);
    return 0.25 * ((lv[index0] + lv[index1]) - (lv[index2] + lv[index3])) / lm[mapd(0, 0, 0, 1)];
  }

  inline RealType dudz(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudz in the cell center by a central difference
    const int index0 = mapd(0, 0, 1, 0);
    const int index1 = mapd(-1, 0, 1, 0);
    const int index2 = mapd(0, 0, -1, 0);
    const int index3 = mapd(-1, 0, -1, 0);
    return 0.25 * ((lv[index0] + lv[index1]) - (lv[index2] + lv[index3])) / lm[mapd(0, 0, 0, 2)];
  }

  inline RealType dvdx(const RealType* const lv, const RealType* const lm) {
    // Evaluate dvdx in the cell center by a central difference
    const int index0 = mapd(1, 0, 0, 1);
    const int index1 = mapd(1, -1, 0, 1);
    const int index2 = mapd(-1, 0, 0, 1);
    const int index3 = mapd(-1, -1, 0, 1);
    return 0.25 * ((lv[index0] + lv[index1]) - (lv[index2] + lv[index3])) / lm[mapd(0, 0, 0, 0)];
  }

  inline RealType dvdz(const RealType* const lv, const RealType* const lm) {
    // Evaluate dvdz in the cell center by a central difference
    const int index0 = mapd(0, 0, 1, 1);
    const int index1 = mapd(0, -1, 1, 1);
    const int index2 = mapd(0, 0, -1, 1);
    const int index3 = mapd(0, -1, -1, 1);
    return 0.25 * ((lv[index0] + lv[index1]) - (lv[index2] + lv[index3])) / lm[mapd(0, 0, 0, 2)];
  }

  inline RealType dwdx(const RealType* const lv, const RealType* const lm) {
    // Evaluate dwdx in the cell center by a central difference
    const int index0 = mapd(1, 0, 0, 2);
    const int index1 = mapd(1, 0, -1, 2);
    const int index2 = mapd(-1, 0, 0, 2);
    const int index3 = mapd(-1, 0, -1, 2);
    return 0.25 * ((lv[index0] + lv[index1]) - (lv[index2] + lv[index3])) / lm[mapd(0, 0, 0, 0)];
  }

  inline RealType dwdy(const RealType* const lv, const RealType* const lm) {
    // Evaluate dwdy in the cell center by a central difference
    const int index0 = mapd(0, 1, 0, 2);
    const int index1 = mapd(0, 1, -1, 2);
    const int index2 = mapd(0, -1, 0, 2);
    const int index3 = mapd(0, -1, -1, 2);
    return 0.25 * ((lv[index0] + lv[index1]) - (lv[index2] + lv[index3])) / lm[mapd(0, 0, 0, 1)];
  }

  inline RealType computeSQuad2D(const RealType* const localVelocity, const RealType* const localMeshsize) {
    return 0.25
           * (4 * dudx(localVelocity, localMeshsize) * dudx(localVelocity, localMeshsize) 
           + 2.0 * (dudy(localVelocity, localMeshsize) + dvdx(localVelocity, localMeshsize)) * (dudy(localVelocity, localMeshsize) + dvdx(localVelocity, localMeshsize)) 
           + 4.0 * dvdy(localVelocity, localMeshsize) * dvdy(localVelocity, localMeshsize));
  }

  inline RealType computeSQuad3D(const RealType* const localVelocity, const RealType* const localMeshsize) {
    return 0.25
           * (4 * dudx(localVelocity, localMeshsize) * dudx(localVelocity, localMeshsize) 
           + 4.0 * dvdy(localVelocity, localMeshsize) * dvdy(localVelocity, localMeshsize) 
           + 4.0 * dwdz(localVelocity, localMeshsize) * dwdz(localVelocity, localMeshsize) 
           + 2.0 * (dudy(localVelocity, localMeshsize) + dvdx(localVelocity, localMeshsize)) * (dudy(localVelocity, localMeshsize) + dvdx(localVelocity, localMeshsize)) 
           + 2.0 * (dudz(localVelocity, localMeshsize) + dwdx(localVelocity, localMeshsize)) * (dudz(localVelocity, localMeshsize) + dwdx(localVelocity, localMeshsize)) 
           + 2.0 * (dvdz(localVelocity, localMeshsize) + dwdy(localVelocity, localMeshsize)) * (dvdz(localVelocity, localMeshsize) + dwdy(localVelocity, localMeshsize)));
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  inline RealType computeF2D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const RealType* const localViscosity, const Parameters& parameters, RealType dt
  ) {
    const RealType Re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (((2.0)*F1(localViscosity, localMeshsize, localVelocity, Re) + F2(localViscosity, localMeshsize, localVelocity, Re)) 
            - du2dx(localVelocity, parameters, localMeshsize) 
            - duvdy(localVelocity, parameters, localMeshsize) 
            + parameters.environment.gx);
  }


  inline RealType computeG2D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const RealType* const localViscosity, const Parameters& parameters, RealType dt
  ) {
    const RealType Re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * ((G1(localViscosity, localMeshsize, localVelocity, Re) + (2.0)*G2(localViscosity, localMeshsize, localVelocity, Re)) 
            - duvdx(localVelocity, parameters, localMeshsize)
            - dv2dy(localVelocity, parameters, localMeshsize) 
            + parameters.environment.gy);
  }

  // TODO check 3D functions
  inline RealType computeF3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const RealType* const localViscosity, const Parameters& parameters, RealType dt
  ) {

    const RealType Re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * ((2.0*F1(localViscosity, localMeshsize, localVelocity, parameters.flow.Re) 
              + F2(localViscosity, localMeshsize, localVelocity, Re) 
              + F3(localViscosity, localMeshsize, localVelocity, Re) )
            - du2dx(localVelocity, parameters, localMeshsize) 
            - duvdy(localVelocity, parameters, localMeshsize)
            - duwdz(localVelocity, parameters, localMeshsize) 
            + parameters.environment.gx);
  }

  inline RealType computeG3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const RealType* const localViscosity, const Parameters& parameters, RealType dt
  ) {
    const RealType Re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * ((G1(localViscosity, localMeshsize, localVelocity, parameters.flow.Re) 
              + 2.0*G2(localViscosity, localMeshsize, localVelocity, Re) 
              + G3(localViscosity, localMeshsize, localVelocity, Re) )
            - dv2dy(localVelocity, parameters, localMeshsize) 
            - duvdx(localVelocity, parameters, localMeshsize)
            - dvwdz(localVelocity, parameters, localMeshsize) 
            + parameters.environment.gx);
  }

  inline RealType computeH3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const RealType* const localViscosity, const Parameters& parameters, RealType dt
  ) {
    const RealType Re = parameters.flow.Re;
    return localVelocity[mapd(0, 0, 0, 2)]
        + dt * ((H1(localViscosity, localMeshsize, localVelocity, Re) 
              + H2(localViscosity, localMeshsize, localVelocity, Re) 
              + 2.0*H3(localViscosity, localMeshsize, localVelocity, Re) )
            - dw2dz(localVelocity, parameters, localMeshsize) 
            - duwdx(localVelocity, parameters, localMeshsize)
            - dvwdy(localVelocity, parameters, localMeshsize) + parameters.environment.gz);
  }


} // namespace Stencils
