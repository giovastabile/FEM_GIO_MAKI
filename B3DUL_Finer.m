function Finer = B3DUL_Finer(Mat,x0,S2,d,ddot,ddot2)
%   Finert = B3DuL_Finert([myy,J1,J2,J3],x0,S2,d,ddot,ddot2)
%        
%   Laskee Reissnerin 3D-palkin hitausvoimavektorin käyttäen 2-pisteen Gaussin integrointia. 
%   Palkissa käytetään lineaarista interpolointia.  Geradin&Cardona-kirjan kaava (6.88)         
%
%
% INPUT:
%   myy     pituustiheys [kg/m], skalaari       
%   J1      hitausmomentti x-akselin ympäri, skalaari [kgm^2/m = kgm]
%   J2      hitausmomentti y-akselin ympäri, skalaari [kgm]
%   J3      hitausmomentti z-akselin ympäri, skalaari (kgm]
%   x0      solmujen alkuasemavektori, dim(x0)=6, x0=[x(1),y(1), z(1), x(2),y(2), z(2)]
%   S2      paikallisen y-akselin suuntavektori, ei tarvitse olla ortogonaalinen palkin x-suunnan kanssa 
%   d       solmusiirtymävektori, dim(d)=12, d=[u1,v1,w1,fii_x1,fii_y1,fii_z1, u2,v2,w2,fii_x2,fii_y2,fii_z2] 
%   ddot    solmunopeusvektori, dim(ddot)=12
%   ddot2   solmukiihtyvyysvektori, dim(ddot2)=12
%
% OUTPUT:
%   Finert  hitausvoimavektori Finer = M*qdot2 + Fgyr,  dim(Finert)=12x1
%

d = d(:); ddot = ddot(:); ddot2 = ddot2(:); 
S2 = S2(:);
x0 = x0(:);
I3 = eye(3,3);
Z3 = zeros(3,3);

myy = Mat(1); J1 = Mat(2); J2 = Mat(3); J3 = Mat(4);

% palkin pituuden laskenta
L = norm(x0(4:6)-x0(1:3));

% palkin S2inen yksikkövektori
E1 = 1/L*(-x0(1:3) + x0(4:6));

% kimmomatriisi R0ER0t (R^6x6)
E2 = S2 - E1'*S2*E1;
E2 = E2/norm(E2);
E3 = cross(E1,E2);
R0 = [E1,E2,E3];
R0JR0t = R0*diag([J1,J2,J3])*R0';

% Gaussin pisteiden koordinaatit
k1= (1-1/sqrt(3))/2;  % = 0.2113
k2= (1+1/sqrt(3))/2;  % = 0.7887

% kiertymävektori ja sen aikaderivaatat Gaussin pisteissä
Y1 = k2*d(4:6) + k1*d(10:12);               % 1. Gaussin pisteessä
Y2 = k1*d(4:6) + k2*d(10:12);               % 2. Gaussin pisteessä
Y1dot = k2*ddot(4:6) + k1*ddot(10:12);      % 1. Gaussin pisteessä
Y2dot = k1*ddot(4:6) + k2*ddot(10:12);      % 2. Gaussin pisteessä
Y1dot2 = k2*ddot2(4:6) + k1*ddot2(10:12);   % 1. Gaussin pisteessä
Y2dot2 = k1*ddot2(4:6) + k2*ddot2(10:12);   % 2. Gaussin pisteessä

% tangenttioperaattorin T arvot Gaussin pisteissä
T1 =  MB_T(Y1);                             % 1. Gaussin pisteessä
T2 =  MB_T(Y2);                             % 2. Gaussin pisteessä
T1dot =  MB_Tdot(Y1dot, Y1);                % 1. Gaussin pisteessä
T2dot =  MB_Tdot(Y2dot, Y2);                   % 2. Gaussin pisteessä

% translaatiokiihtyvyysvektorin arvot Gaussin pisteissä
x1dot2 = k2*ddot2(1:3) + k1*ddot2(7:9);     % 1. Gaussin pisteessä
x2dot2 = k1*ddot2(1:3) + k2*ddot2(7:9);     % 2. Gaussin pisteessä

% aineellisen kulmanopeuden ja sen aikaderivaatan arvot Gaussin pisteissä
Omega1 = T1*Y1dot;                          % 1. Gaussin pisteessä
Omega2 = T2*Y2dot;                          % 2. Gaussin pisteessä
Omega1dot = T1*Y1dot2 + T1dot*Y1dot;        % 1. Gaussin pisteessä
Omega2dot = T2*Y2dot2 + T2dot*Y2dot;        % 2. Gaussin pisteessä

% interpolaatiomatriisin P arvot gaussin pisteissä 
P1 = [k2*I3,Z3,k1*I3,Z3; Z3,k2*I3,Z3,k1*I3];   % 1. Gaussin pisteessä
P2 = [k1*I3,Z3,k2*I3,Z3; Z3,k1*I3,Z3,k2*I3];   % 2. Gaussin pisteessä

Finer = L/2*( P1'*[myy*x1dot2;  T1'*(R0JR0t*Omega1dot + MB_mato(Omega1)*R0JR0t*Omega1)] + P2'*[myy*x2dot2;  T2'*(R0JR0t*Omega2dot + MB_mato(Omega2)*R0JR0t*Omega2)]);