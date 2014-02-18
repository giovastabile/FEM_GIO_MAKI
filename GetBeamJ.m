function J=GetBeamJ(Roo,Iy,Iz)
% ****************************************************************************************************************
% It calculates the inertia tensor of the beam
%       INPUT:
% 
%       Iy          Moment of inertia respect to y         
%       Iz          Moment of inertia respect to z  
%       Roo         Scalar Mass Density
%
%       OUTPUT:
% 
%       J           Inertia tensor
%       
% 
%       J=Roo*[Jt  0   0]
%             [0   Iy  0]
%             [0   0  Iz]    
% ****************************************************************************************************************
  J(1,1)=Roo*(Iy+Iz);
  J(2,2)=Roo*Iy;
  J(3,3)=Roo*Iz;