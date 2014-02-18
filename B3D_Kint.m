function [Kmat,Kgeom] = B3D_Kint(EA,GA2,GA3,GJ,EI2,EI3,x0,suunta,d)
%   [Kmat,Kgeom] = B3D_Kint(EA,GA2,GA3,GJ,EI2,EI3,x0,suunta,d)
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
%
% OUTPUT:
%   Kmat     materiaalinen j�ykkyysmatriisi, dim(Kmat)=12x12
%   Kgeom    geometrinen j�ykkyysmatriisi, dim(Kgeom)=12x12
%	

      
%
%   tarkistettu 3.11.2001, Jari M�kinen, <jari@mohr.me.tut.fi>
%

d = d(:);
x0 = x0(:);
suunta = suunta(:);
I3 = eye(3,3);
Z3 = zeros(3,3);

% 1. solmun kiertym�vektori Y1
Y1 = d(4:6);
% 2. solmun kiertym�vektori Y2 
Y2 = d(10:12); 
% 2. solmun kiertym�kulma fii_Y2
fii_Y2 = norm(Y2);

% vaihto on totuusarvoinen muuttuja 
vaihto = 0;
if fii_Y2 >  0.01,  
    Y2c = Y2-2*pi/fii_Y2*Y2; 
    t13 =  norm(Y2-Y1) - norm(Y2c-Y1);
    % mik�li t13>0, niin 2. solmuun kiertym�vektorin kutistus
    if t13  > 0,  
        vaihto = 1;  
        % kiertym�vektorin vaihtoa vastaava kinemaattimen matriisi Bvaihto e R^3x3
        Bvaihto = (1-2*pi/fii_Y2)*I3+2*pi/fii_Y2^3*(Y2*Y2');
        % p�ivitet��n Y2-vektori
        Y2 = Y2c; 
    else
        vaihto = 0;
    end
end
  
% pieni kulma: kun kulma on t�t� pienempi lasketaan sis�inen voimavektori yksinkertaisemman kaavan mukaisesti
fii_epsi = 1e-6; 
% palkin pituuden laskenta
L = norm(x0(4:6)-x0(1:3));

% kiertym�suureet gaussin pisteess�
Yg = 1/2*(Y1 + Y2);
% kiertym�suureiden derivaatat gaussin pisteess�
Yg_s = 1/L*(Y2 - Y1); 

% rotaatiomatriisi R, transformaatiomatriisit T 
fii = norm(Yg);  % kulma fii
R = MB_R(Yg);
T = MB_T(Yg);

% kimmoviivan derivaatta gaussin pisteess� 
x0dp_dot = 1/L*(x0(4:6)+d(7:9)-x0(1:3)-d(1:3)); 

% vektori transpoosi(R)*xdot
Rtxdot = R'*x0dp_dot;

% palkin suuntainen yksikk�vektori
E1 = 1/L*(-x0(1:3) + x0(4:6));

% venym�vektori GAMMA ;
G = Rtxdot - E1;

% Kaareutumavektori K (matriisi T kertaa vektori Yg_s
K = T*Yg_s;

% kimmomatriisi R0ER0t (R^6x6)
E2 = suunta - E1'*suunta*E1;
E2 = E2/norm(E2);
E3 = cross(E1,E2);
R0 = [E1,E2,E3];
R0ER0t = [R0*diag([EA,GA2,GA3])*R0', zeros(3,3); zeros(3,3), R0*diag([GJ,EI2,EI3])*R0'];

% sis�inen normaalivoima ja momentti NjaM 
N = R0*diag([EA,GA2,GA3])*R0'*G;
M = R0*diag([GJ,EI2,EI3])*R0'*K;

C1 = MB_C1(Yg_s,Yg);
Bh = [R', Z3, mato(Rtxdot)*T; Z3, T, C1];

Lambda = [-1/L*eye(6,6), 1/L*eye(6,6); Z3, 1/2*I3, Z3, 1/2*I3];
B = Bh*Lambda;
C2_1 = MB_C2(M,Yg);
C2_2 = MB_C2(mato(N)*Rtxdot,Yg);
C3 = MB_C3(M,Yg_s,Yg); 

% materiaalinen j�ykkyysmatriisi
Kmat = L*B'*R0ER0t*B;

% geometrinen j�ykkyyematriisi
Kg = [Z3,Z3,-R*mato(N)*T; Z3,Z3,C2_1; T'*mato(N)*R', C2_1', C3+C2_2+T'*mato(N)*mato(Rtxdot)*T]; 
Kgeom = L*Lambda'*Kg*Lambda;

if vaihto,
    % sis�inen voimavektori 
    F_int = Beam3D_Fint(EA,GA2,GA3,GJ,EI2,EI3,x0,suunta,[d(1:9); Y2]);
    Y2 = d(10:12); 
    
    % geometrinen j�ykkyysmatriisi
    BB = [eye(9,9), zeros(9,3); zeros(3,9), Bvaihto];
    Kg = BB'*Kgeom*BB;
    e = Y2/norm(Y2);
    Kg1= 2*pi/fii_Y2^2*( F_int(10:12)*e' + e*F_int(10:12)' - 3*(F_int(10:12)'*e)*(e*e') + (e'*F_int(10:12))*I3);
    Kgeom = Kg + [zeros(9,12); [zeros(3,9), Kg1]];
    
    %materiaalinen j�ykkyysmatriisi    
    Kmat = BB'*Kmat*BB;
end    
