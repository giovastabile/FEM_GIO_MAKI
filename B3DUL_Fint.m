function [Fint,R,K] = B3DUL_Fint(Mat,x0,suunta,d,Rref, Kref)
%   [Fint,R,K] = B3DuL_Fint([EA,GA2,GA3,GJ,EI2,EI3],x0,suunta,d,Rref, Kref)
% 
%   It Calculates the Simo-Reiossner 3d internal forces of the beam using a
%   1-point Gauss integration rule with linear interpolation function.
%
% INPUT:
%   Mat     [EA GA2 GA3 GJ EI2 EI3] 
% 
%   with:
% 
%   EA      Axial Stiffness    
%   GA2     Shear Stiffness along y direction
%   GA3     Shear Stiffness along z direction
%   GJ      Torsional Stiffness (Rotation respect to x)
%   EI2     Bending Stiffness   (Rotation respect to y)
%   EI3     Bending Stiffness   (Rotation respect to y)
% 
%   x0      Element nodes in the initial position, dim(x0)=6, x0=[x(1),y(1), z(1), x(2),y(2), z(2)]
%   suunta  y axis of the local direction. It doesn't need to be orthogonal to the x-axis 
%   d       solmusiirtymävektori, dim(d)=12, d=[u(1),v(1), w(1), fii_x(1), fii_y(1), fii_z(1), u(2),v(2), w(2), fii_x(2), fii_y(2), fii_z(2)] 
%   Rref    Reference system, dim(Rref)=3x3
%   Kref    State of curvature, dim(Kref)=3x1
%
% OUTPUT:
%   Fint    Internal force vector, dim(Fint)=12x1
%


EA=Mat(1);
GA2=Mat(2);
GA3=Mat(3);
GJ=Mat(4);
EI2=Mat(5);
EI3=Mat(6);


d = d(:);
x0 = x0(:);
suunta = suunta(:);
I3 = eye(3,3);
Z3 = zeros(3,3);


% 1. solmun kiertymävektori Y1
Y1 = d(4:6);
% 2. solmun kiertymävektori Y2 
Y2 = d(10:12); 
  
% alkutilassa palkin pituus 
L = norm(x0(4:6)-x0(1:3));


% kiertymävektori Gaussin pisteessä
Yg = 1/2*(Y1 + Y2);
% kiertymävektorin paikkaderivaata Gaussin pisteessä
Yg_s = 1/L*(Y2 - Y1); 


% Gaussin pisteen kiertymämatriisi R=Rref*Rinc
Rinc = MB_R(Yg);
R = Rref*Rinc;


% transformaatiomatriisi T 
T = MB_T(Yg);


% kimmoviivan derivaatta Gaussin pisteessä 
x0dp_dot = 1/L*(x0(4:6)+d(7:9)-x0(1:3)-d(1:3)); 


% vektori transpoosi(R)*xdot
Rtxdot = R'*x0dp_dot;


% palkin suuntainen yksikkövektori
E1 = 1/L*(-x0(1:3) + x0(4:6));


% venymävektori GAMMA, Gaussin pisteessä
G = Rtxdot - E1;


% Gaussin pisteen kaarevuusvektori  K = Rinc'*Kref + Kinc
Kinc = T*Yg_s;
K = Rinc'*Kref + Kinc;


% alkutilan kiertymämatriisi R0 (R^6x6)
E2 = suunta - E1'*suunta*E1;  % ortonormeerataan E2
E2 = E2/norm(E2);
E3 = cross(E1,E2);
R0 = [E1,E2,E3];


% sisäinen normaalivoima ja momentti NjaM Gaussin pisteessä 
N = R0*diag([EA,GA2,GA3])*R0'*G;
M = R0*diag([GJ,EI2,EI3])*R0'*K;


% T', T:n paikkaderivaatta Gaussin pisteessä
Tprime= MB_Tdot(Yg_s,Yg);


% interpolaatiomatriisi Q, (6,77) Gaussin pisteessä
Q = [-1/L*eye(6,6), 1/L*eye(6,6); Z3, 1/2*I3, Z3, 1/2*I3];
% kinemaattinen matriisi, (6.79) 
B = [R', Z3, MB_mato(Rtxdot)*T; Z3, T, MB_mato(K)*T + Tprime]*Q;


% sisäinen solmuvoimavektori, lasketaan yhdellä Gaussin pisteellä, (6.81)
Fint = L*B'*[N;M]; 