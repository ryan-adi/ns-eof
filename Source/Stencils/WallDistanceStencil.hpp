#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil to compute the velocity once the pressure has been found.
   */
  class WallDistanceStencil: public FieldStencil<TurbulentFlowField> {
  private: 
  const RealType X_max_;
  const RealType Y_max_;
  const RealType Z_max_;
  const RealType x_step_;
  const RealType y_step_;
  
  public:
    WallDistanceStencil(const Parameters& parameters);
    ~WallDistanceStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) ; 
    void apply(TurbulentFlowField& flowField, int i, int j, int k) ; 
  };

} // namespace Stencils
