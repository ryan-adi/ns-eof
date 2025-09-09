#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil to compute the turbulent viscosity.
   */
  class TVStencil: public FieldStencil<TurbulentFlowField> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    TVStencil(const Parameters& parameters);
    ~TVStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j);
    void apply(TurbulentFlowField& flowField, int i, int j, int k);
  };

} // namespace Stencils
