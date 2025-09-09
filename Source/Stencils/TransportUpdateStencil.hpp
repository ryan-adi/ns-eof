#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil to compute the turbulent vorticty omega
   */
  class TransportUpdateStencil: public FieldStencil<TurbulentFlowField> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.

  public:
    TransportUpdateStencil(const Parameters& parameters);
    ~TransportUpdateStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j);
    void apply(TurbulentFlowField& flowField, int i, int j, int k);
  };

} // namespace Stencils
