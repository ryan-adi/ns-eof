#include "StdAfx.hpp"

#include "MaxTVStencil.hpp"

Stencils::MaxTVStencil::MaxTVStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  BoundaryStencil<TurbulentFlowField>(parameters) {

  reset();
}

void Stencils::MaxTVStencil::apply(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxTVStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxTVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxTVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxTVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }


void Stencils::MaxTVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTVStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTVStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxTVStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j) {
  if (fabs(flowField.getTurbViscosity().getScalar(i, j)) > maxValues_[0]) {
    maxValues_[0] = fabs(flowField.getTurbViscosity().getScalar(i, j));
  }
}

void Stencils::MaxTVStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j, int k) {
  if (fabs(flowField.getTurbViscosity().getScalar(i, j)) > maxValues_[0]) {
    maxValues_[0] = fabs(flowField.getTurbViscosity().getScalar(i, j, k));
  }
}
void Stencils::MaxTVStencil::reset() { maxValues_[0] = 0; }

const RealType* Stencils::MaxTVStencil::getMaxValues() const { return maxValues_; }
