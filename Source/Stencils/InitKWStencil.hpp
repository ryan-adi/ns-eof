#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil to compute the inital values for k and omega
   */
  class InitKWStencil: public FieldStencil<TurbulentFlowField> {
  private:
    const RealType xLimit_; //! Size of step in x-direction
    const RealType yLimit_; //! Same as for x
  
  public:
    InitKWStencil(const Parameters& parameters);
    ~InitKWStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) ; 
    void apply(TurbulentFlowField& flowField, int i, int j, int k) ; 
  };

} // namespace Stencils
