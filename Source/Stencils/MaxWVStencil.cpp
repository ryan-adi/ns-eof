#include "StdAfx.hpp"

#include "MaxWVStencil.hpp"

Stencils::MaxWVStencil::MaxWVStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  BoundaryStencil<TurbulentFlowField>(parameters) {

  reset();
}

void Stencils::MaxWVStencil::apply(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxWVStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxWVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxWVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxWVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxWVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }


void Stencils::MaxWVStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxWVStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxWVStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxWVStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxWVStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxWVStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxWVStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j) {
  if (fabs(flowField.getTurbOmega().getScalar(i, j)) > maxValues_[0]) {
    maxValues_[0] = fabs(flowField.getTurbOmega().getScalar(i, j));
  }
}

void Stencils::MaxWVStencil::cellMaxValue(TurbulentFlowField& flowField, int i, int j, int k) {
  if (fabs(flowField.getTurbOmega().getScalar(i, j)) > maxValues_[0]) {
    maxValues_[0] = fabs(flowField.getTurbOmega().getScalar(i, j, k));
  }
}
void Stencils::MaxWVStencil::reset() { maxValues_[0] = 0; }

const RealType* Stencils::MaxWVStencil::getMaxValues() const { return maxValues_; }