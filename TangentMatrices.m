function [M,C,K,B] = TangentMatrices(t,y0,v0,a0,l0)
% funktio palauttaa tehtävän tangenttimatriisit

  C=Damping(t,y0,v0,a0); 
  M=Mass(t,y0,v0,a0);  
  K=Stiffness_UL(t,y0,v0,a0,l0);  
  B=ConstraintG(t,y0,v0,a0);
  
  % end function
