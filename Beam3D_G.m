function Fg=Beam3D_G(X, rho)
%***************************************************************
%Beam element gravity load
%
%Parameters     X   Element initial position
%               rho Element linear density
%
%Returns        Fg() load vector 12 * 1
%
%Programmer:    Heikki Marjamäki
%Date:          1.10.2004
%
%Modified:      
%
%***************************************************************
%  disp('Beam3D_G')

global FEM;

  C_HALF=1/2;
  x21 = X(4)-X(1);
  y21 = X(5)-X(2);
  z21 = X(6)-X(3);
  l02 = x21 * x21 + y21*y21+z21*z21;
  L0=sqrt(l02);
  
  
  Ele_m = rho * L0;
  %
  Fg(1, 1) = C_HALF * Ele_m*FEM.PARA.Gx;
  Fg(2, 1) = C_HALF * Ele_m*FEM.PARA.Gy;
  Fg(3, 1) = C_HALF * Ele_m*FEM.PARA.Gz;
  Fg(7, 1) = C_HALF * Ele_m*FEM.PARA.Gx;
  Fg(8, 1) = C_HALF * Ele_m*FEM.PARA.Gy;
  Fg(9, 1) = C_HALF * Ele_m*FEM.PARA.Gz;
  Fg(12, 1) =0;
  % end Function
