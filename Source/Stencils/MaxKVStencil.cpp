#include "StdAfx.hpp"

#include "MaxKVStencil.hpp"

Stencils::MaxKVStencil::MaxKVStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  BoundaryStencil<TurbulentFlowField>(parameters) {

  reset();
}

void Stencils::MaxKVStencil::apply(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxKVStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxKVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxKVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxKVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxKVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }


void Stencils::MaxKVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxKVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxKVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxKVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxKVStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxKVStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxKVStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j) {
  if (fabs(flowField.getTurbKineticEnergy().getScalar(i, j)) > maxValues_[0]) {
    maxValues_[0] = fabs(flowField.getTurbKineticEnergy().getScalar(i, j));
  }
}

void Stencils::MaxKVStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j, int k) {
  if (fabs(flowField.getTurbKineticEnergy().getScalar(i, j)) > maxValues_[0]) {
    maxValues_[0] = fabs(flowField.getTurbKineticEnergy().getScalar(i, j, k));
  }
}
void Stencils::MaxKVStencil::reset() { maxValues_[0] = 0; }

const RealType* Stencils::MaxKVStencil::getMaxValues() const { return maxValues_; }
