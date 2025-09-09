#pragma once

#include "BoundaryStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** A stencil to set the input velocity in channel flows. Only implements the applyLeftWall(...) methods.
   */
  class BFInputVelocityStencil: public BoundaryStencil<FlowField> {
  private:
    RealType stepSize_; //! Fixes the size of the step. If zero, is just channel flow.

  public:
    BFInputVelocityStencil(const Parameters& parameters);
    ~BFInputVelocityStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;
  };

  /** FGH stencil which corresponds to one-sided Dirichlet conditions (i.e. the channel flow profile).
   *  Only implements the applyLeftWall(...) methods.
   */
  class BFInputFGHStencil: public BoundaryStencil<FlowField> {
  private:
    const RealType stepSize_; //! Fixes the size of the step. If zero, is just channel flow.

  public:
    BFInputFGHStencil(const Parameters& parameters);
    ~BFInputFGHStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;
  };

  /** TV stencil which corresponds to one-sided Dirichlet conditions (i.e. the channel flow profile).
   *  Only implements the applyLeftWall(...) methods.
   */
  class BFInputTVStencil: public BoundaryStencil<TurbulentFlowField> {
  private:
    const RealType stepSize_; //! Fixes the size of the step. If zero, is just channel flow.

  public:
    BFInputTVStencil(const Parameters& parameters);
    ~BFInputTVStencil() override = default;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

  /** K stencil which corresponds to one-sided Dirichlet conditions (i.e. the channel flow profile).
   *  Only implements the applyLeftWall(...) methods.
   */
  class BFInputKStencil: public BoundaryStencil<TurbulentFlowField> {
  private:
    const RealType stepSize_; //! Fixes the size of the step. If zero, is just channel flow.

  public:
    BFInputKStencil(const Parameters& parameters);
    ~BFInputKStencil() override = default;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

  /** W stencil which corresponds to one-sided Dirichlet conditions (i.e. the channel flow profile).
   *  Only implements the applyLeftWall(...) methods.
   */
  class BFInputWStencil: public BoundaryStencil<TurbulentFlowField> {
  private:
    const RealType stepSize_; //! Fixes the size of the step. If zero, is just channel flow.

  public:
    BFInputWStencil(const Parameters& parameters);
    ~BFInputWStencil() override = default;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
