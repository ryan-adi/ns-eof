#include "StdAfx.hpp"

#include "WallDistanceStencil.hpp"

/**
 * WallDistanceStencil 
 * 
*/

Stencils::WallDistanceStencil::WallDistanceStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters), 
    X_max_(parameters.geometry.lengthX),
    Y_max_(parameters.geometry.lengthY),
    Z_max_(parameters.geometry.lengthZ),
    x_step_(parameters.bfStep.xRatio * parameters.geometry.lengthX),
    y_step_(parameters.bfStep.yRatio * parameters.geometry.lengthY){}

void Stencils::WallDistanceStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const int      obstacle = flowField.getFlags().getValue(i, j);
  // Load local velocities into the center layer of the local array RA: TODO use local velocity ???
  
  const RealType posX     = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i,j); // x-pos left bottom front node
  const RealType posY     = parameters_.meshsize->getPosY(i, j) + 0.5 * parameters_.meshsize->getDy(i,j); // y-pos left bottom front node
  //const RealType lX       = parameters_.geometry.lengthX;
  ScalarField& wallDistance  = flowField.getNearestWallH();   



if ((obstacle & OBSTACLE_SELF) == 0) {
  //apply cavity method
  if(parameters_.simulation.scenario == "cavity"){
    //check what is the minimum distance between the minimum distance from the wall 
    //in both orientations for both directions
     wallDistance.getScalar(i, j) =std::min(std::min((Y_max_-posY),posY),std::min((X_max_-posX),posX));


  //apply channel method    
  }else if(parameters_.simulation.scenario == "channel" && x_step_ == 0){
    //check y direction minimum distance from wall
    wallDistance.getScalar(i, j) = std::min((Y_max_ - posY),posY);


  //apply backwardfacing step method  
  }else if(parameters_.simulation.scenario == "channel" && x_step_ != 0){

    //if cell is before the step than only y direction is considered
    //and the geometry becomes channel with height reduced by y_step_
    if(posX < x_step_){
      wallDistance.getScalar(i, j) =  std::min((posY - y_step_),(Y_max_ - posY));
      //if cell is after the step than also x direction and diagonal distance
      //must be considered
    }else if(posX >= x_step_){
      //checking if it's over or below the step height
      if (posY <= y_step_){
        // get minimum between horizontal distance from the step and 
        //distance from the the ground
        wallDistance.getScalar(i, j) =  std::min((posY),std::abs(posX - x_step_));

        //if point is above the step near the ending part of the channel
        //must check also diagonal distance other than vertical
        //thanks to last if no x direction has to be checked
      }else if (posY > y_step_){
      //compute diagonal distance from the corner
      //by construction corner point has coordinates(x_step_,y_step_y)
      RealType diag = std::sqrt(std::pow(posX-x_step_,2)+std::pow(posY-y_step_,2));
      //take minimum distance between :
      //distance from ceiling->(Y_max_-posY)
      //diagonal distance from the corner point->diag
      //distance from the floor->PosY
      wallDistance.getScalar(i, j) =  std::min(std::min(posY,(Y_max_-posY)), diag);
    } 
  
    }

}
}else{
   wallDistance.getScalar(i, j) = 0.0;//MY_FLOAT_MAX;
}

}
// TODO 3D case not done
void Stencils::WallDistanceStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  const int      obstacle = flowField.getFlags().getValue(i, j,k);
  // Load local velocities into the center layer of the local array RA: TODO use local velocity ???
  
  const RealType posX     = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i,j,k); // x-pos left bottom front node
  const RealType posY     = parameters_.meshsize->getPosY(i, j, k) + 0.5 * parameters_.meshsize->getDy(i,j,k); // y-pos left bottom front node
  const RealType posZ     = parameters_.meshsize->getPosZ(i, j, k) + 0.5 * parameters_.meshsize->getDz(i,j,k); // z-pos left bottom front node
  //const RealType lX       = parameters_.geometry.lengthX;
  ScalarField& wallDistance  = flowField.getNearestWallH();   

  RealType dmin = 0;

if ((obstacle & OBSTACLE_SELF) == 0) {
  //apply cavity method
  if(parameters_.simulation.scenario == "cavity"){
    //check what is the minimum distance between the minimum distance from the wall 
    //in both orientations for both directions
     dmin = std::min(std::min((Y_max_-posY),posY),std::min((X_max_-posX),posX));
     wallDistance.getScalar(i, j, k) =std::min(dmin,std::min((Z_max_-posZ),posZ));
  //apply channel method    
  }else if(parameters_.simulation.scenario == "channel" && x_step_ == 0){
    //check y direction minimum distance from wall

     dmin = std::min((Y_max_ - posY),posY);
     wallDistance.getScalar(i, j, k) = std::min(dmin,std::min((Z_max_ - posZ),posZ));

  //apply backwardfacing step method  
  }else if(parameters_.simulation.scenario == "channel" && x_step_ != 0){

    //if cell is before the step than only y direction is considered
    //and the geometry becomes channel with height reduced by y_step_
    if(posX < x_step_){
      wallDistance.getScalar(i, j, k) =  std::min(std::min((posY - y_step_),(Y_max_ - posY)),std::min((Z_max_ - posZ),posZ));
      //if cell is after the step than also x direction and diagonal distance
      //must be considered
    }else if(posX >= x_step_){
      //checking if it's over or below the step height
      if (posY <= y_step_){
        // get minimum between horizontal distance from the step and 
        //distance from the the ground
        wallDistance.getScalar(i, j, k) =  std::min(std::min((posY),std::abs(posX - x_step_)),std::min((Z_max_ - posZ),posZ));
        //if point is above the step near the ending part of the channel
        //must check also diagonal distance other than vertical
        //thanks to last if no x direction has to be checked
      }else if (posY > y_step_){
      //compute diagonal distance from the corner
      //by construction corner point has coordinates(x_step_,y_step_y)
      RealType diag = std::sqrt(std::pow(posX-x_step_,2)+std::pow(posY-y_step_,2));
      //take minimum distance between :
      //distance from ceiling->(Y_max_-posY)
      //diagonal distance from the corner point->diag
      //distance from the floor->PosY
      wallDistance.getScalar(i, j, k) = std::min(std::min(std::min(posY,(Y_max_-posY)), diag),std::min((Z_max_ - posZ),posZ));
    } 
  
    }

}
}else{
   wallDistance.getScalar(i, j, k) = 0.0;//MY_FLOAT_MAX;
}
 
};




