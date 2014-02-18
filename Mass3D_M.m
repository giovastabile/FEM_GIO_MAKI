function M=Mass2D_M(EleMass)
%***************************************************************
%Point mass element mass matrix
%
%Parameters     EleMass  Element mass
%
%Returns        M() mass matrix 6 * 6
%
%Programmer:    Heikki Marjamäki
%Date:          19.4.2005
%
%Modified:      
%
%***************************************************************
  %Mass matrix
  M=zeros(6,6);
  M(1:3,1:3) = EleMass*eye(3);
  
end %Function
