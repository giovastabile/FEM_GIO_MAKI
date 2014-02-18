function [Kmat,Kgeom] = B3DUL_Kint(Mat,x0,suunta,d,Rref, Kref)
%   [Kmat,Kgeom] = B3DuL_Kint([EA,GA2,GA3,GJ,EI2,EI3],x0,suunta,d,Rref, Kref)
%        
%   Laskee Reissnerin 3D-palkin sis�isen j�ykkyysmatriisin k�ytt�en 1-pisteen gaussin integrointia, 
%   jossa  lineaariset interpolaatiofuktiot.    
%
% INPUT:
%   EA      veto-puristusj�ykkyys, skalaari     
%   GA2     leikkausj�ykkyys y-suuntaan, skalaari
%   GA3     leikkausj�ykkyys z-suuntaan, skalaari
%   GJ      v��nt�j�ykkyys x-akselin ymp�ri, skalaari
%   EI2     taivutusj�ykkyys y-akselin ymp�ri, skalaari
%   EI3     taivutusj�ykkyys z-akselin ymp�ri, skalaari
%   x0      solmujen alkuasemavektori, dim(x0)=6, x0=[x(1),y(1), z(1), x(2),y(2), z(2)]
%   suunta  paikallisen y-akselin suuntavektori, ei tarvitse olla ortogonaalinen palkin x-suunnan kanssa 
%   d       solmusiirtym�vektori, dim(d)=12, d=[u(1),v(1), w(1), fii_x(1), fii_y(1), fii_z(1), u(2),v(2), w(2), fii_x(2), fii_y(2), fii_z(2)] 
%   Rref    vertailutilan kiertym�matriisi gaussin pisteess�, dim(Rref)=3x3
%   Kref    vertailutilan kaarevuusvektori gaussin pisteess�, dim(Kref=3x1
%
% OUTPUT:
%   Kmat     materiaalinen j�ykkyysmatriisi, dim(Kmat)=12x12
%   Kgeom    geometrinen j�ykkyysmatriisi, dim(Kgeom)=12x12
%       

EA=Mat(1);
GA2=Mat(2);
GA3=Mat(3);
GJ=Mat(4);
EI2=Mat(5);
EI3=Mat(6);


% tehd��n vektoreista pystyvektoreita
d = d(:);
x0 = x0(:);
Kref = Kref(:);
suunta = suunta(:);


% apumatriiseja
I3 = eye(3,3);
Z3 = zeros(3,3);


% palkin pituuden laskenta alkutilassa
L = norm(x0(4:6)-x0(1:3));


% kiertym�vektori Gaussin pisteess�
Yg = 1/2*(d(4:6) + d(10:12));
Yg_s = 1/L*(d(10:12) - d(4:6)); 


% Gaussin pisteen kiertym�matriisi R=Rref*Rinc
Rinc = MB_R(Yg);
R = Rref*Rinc;


% kimmoviivan derivaatta Gaussin pisteess� 
x0dp_dot = 1/L*(x0(4:6)+d(7:9)-x0(1:3)-d(1:3)); 


% apuvektori transpoosi(R)*xdot
Rtxdot = R'*x0dp_dot;


% palkin suuntainen yksikk�vektori
E1 = 1/L*(-x0(1:3) + x0(4:6));


% venym�vektori GAMMA ;
G = Rtxdot - E1;


% Gaussin pisteen kaarevuusvektori  K = Rinc'*Kref + Kinc, (Rinc=I, Kinc=0, kun Y->0)
% K = Kref;
Kinc=MB_T(Yg)*Yg_s;
K = Rinc'*Kref + Kinc;

% alkutilan kiertym�matriisi R0 (R^6x6)
E2 = suunta - E1'*suunta*E1;
E2 = E2/norm(E2);
E3 = cross(E1,E2);
R0 = [E1,E2,E3];
R0ER0t = [R0*diag([EA,GA2,GA3])*R0', zeros(3,3); zeros(3,3), R0*diag([GJ,EI2,EI3])*R0'];


% sis�inen normaalivoima ja momentti NjaM Gaussin pisteess�
N = R0*diag([EA,GA2,GA3])*R0'*G;
M = R0*diag([GJ,EI2,EI3])*R0'*K;


% interpolointimatriisi Q, (6.77)
Q = [-1/L*eye(6,6), 1/L*eye(6,6); Z3, 1/2*I3, Z3, 1/2*I3];


% kinemaattinen matriisi B, (6.79)
B = [R', Z3, MB_mato(Rtxdot); Z3, I3, MB_mato(K)]*Q;


% materiaalinen j�ykkyysmatriisi Kmat(0), (6.96)
Kmat = L*B'*R0ER0t*B;


% geometrinen j�ykkyysmatriisi Kgeom(0), (6.121)
W = [Z3,Z3,-R*MB_mato(N); Z3,Z3,-0.5*MB_mato(M); MB_mato(N)*R', 0.5*MB_mato(M), 0.5*(MB_mato(N)*MB_mato(Rtxdot)+MB_mato(Rtxdot)*MB_mato(N)) + 0.5*(MB_mato(M)*MB_mato(K)+MB_mato(K)*MB_mato(M))]; 
Kgeom = L*Q'*W*Q;
