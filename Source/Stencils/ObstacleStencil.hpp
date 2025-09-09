#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Compute all velocities on obstacle cells
   */
  class ObstacleStencil: public FieldStencil<FlowField> {
  public:
    ObstacleStencil(const Parameters& parameters);
    ~ObstacleStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
    
  };

  class TurbObstacleStencil: public FieldStencil<TurbulentFlowField> {
    private:
    RealType computeW2D(RealType Re, RealType Dx);
    public:
      TurbObstacleStencil(const Parameters& parameters);
      ~TurbObstacleStencil() override = default;
      void apply(TurbulentFlowField& flowField, int i, int j) override;
      void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
