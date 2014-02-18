function [FaccA,FaccB,M] = B3D_Facc(Arhoo,J1,J2,J3,x0,suunta,d,ddot,ddot2)
%   [FaccA,FaccB,M] = B3D_Facc(Arhoo,J1,J2,J3,x0,suunta,d,ddot,ddot2)
%	 
%   Laskee Reissnerin 3D-palkin kiihtytyys voimavektorin käyttäen 2-pisteen gaussin integrointia, 
%   jossa on lineaariset interpolaatiofuktiot.  	
%
%
% INPUT:
%   Arhoo   A*rhoo = pituustiheys [kg/m], skalaari	
%   J1      hitausmomentti x-akselin ympäri, skalaari [kgm^2/m = kgm]
%   J2      hitausmomentti y-akselin ympäri, skalaari [kgm]
%   J3      hitausmomentti z-akselin ympäri, skalaari (kgm]
%   x0      solmujen alkuasemavektori, dim(x0)=6, x0=[x(1),y(1), z(1), x(2),y(2), z(2)]
%   suunta  saikallisen y-akselin suuntavektori, ei tarvitse olla ortogonaalinen palkin x-suunnan kanssa 
%   d       solmusiirtymävektori, dim(d)=12, d=[u1,v1,w1,fii_x1,fii_y1,fii_z1, u2,v2,w2,fii_x2,fii_y2,fii_z2] 
%   ddot    solmunopeusvektori, dim(ddot)=12
%   ddot2   solmukiihtyvyysvektori, dim(ddot2)=12
%
% OUTPUT:
%   FaccA   hitausvoimavektori A, M*ddot2, dim(FaccA)=12
%   FaccB   hitausvoimavektori B, dim(FaccB)=12
%   M       massamatriisi dim(M)=12x12
%

%
%   tarkistettu 14.10.2002, Jari Mäkinen, <jari@mohr.me.tut.fi>
%

d = d(:); ddot = ddot(:); ddot2 = ddot2(:); 
suunta = suunta(:);
x0 = x0(:);
I3 = eye(3,3);
Z3 = zeros(3,3);

% 1. solmun kiertymävektori Y1
Y1 = d(4:6);
Y1dot = ddot(4:6);
Y1dot2 = ddot2(4:6);
% 2. solmun kiertymävektori Y2 
Y2 = d(10:12); 
Y2dot = ddot(10:12);
Y2dot2 = ddot2(10:12);
% 2. solmun kiertymäkulma fii_Y2
fii_Y2 = norm(Y2);

% vaihto on totuusarvoinen muuttuja 
vaihto = 0;
if fii_Y2 >  0.01,  
    Y2c = Y2-2*pi/fii_Y2*Y2; 
    t13 =  norm(Y2-Y1) - norm(Y2c-Y1);
    % mikäli t13>0, niin 2. solmuun kiertymävektorin kutistus
    if t13  > 0,  
        vaihto = 1;  
        % kiertymävektorin vaihtoa vastaava kinemaattimen matriisi Bvaihto e R^3x3
        Bvaihto = (1-2*pi/fii_Y2)*I3+2*pi/fii_Y2^3*(Y2*Y2'); % tarvitaan jatkossa
        Bvaihtodot = 2*pi/fii_Y2^3*( (Y2'*Y2dot)*I3 + (Y2dot*Y2' + Y2*Y2dot') - 3/fii_Y2^2*(Y2'*Y2dot)*Y2*Y2' );
        Y2dotc = Bvaihto*Y2dot;
        Y2dot2 = Bvaihto*Y2dot2 + Bvaihtodot*Y2dot;
        % päivitetään Y2-vektori
        Y2 = Y2c; 
        Y2dot = Y2dotc;

    else
        vaihto = 0;
    end
end
  
% pieni kulma: kun kulma on tätä pienempi lasketaan sisäinen voimavektori yksinkertaisemman kaavan mukaisesti
fii_epsi = 1e-6; 
% palkin pituuden laskenta

L = norm(x0(4:6)-x0(1:3));
% palkin suuntainen yksikkövektori
E1 = 1/L*(-x0(1:3) + x0(4:6));

% kiertymäsuureet gaussin pisteessä
Yg = 1/2*(Y1 + Y2);
Ygdot = 1/2*(Y1dot + Y2dot);

% rotaatiomatriisi R, transformaatiomatriisit T 
R = MB_R(Yg);
T = MB_T(Yg);
Tdot = MB_Tdot(Ygdot,Yg);
% kimmomatriisi R0ER0t (R^6x6)
E2 = suunta - E1'*suunta*E1;
E2 = E2/norm(E2);
E3 = cross(E1,E2);
R0 = [E1,E2,E3];
R0JR0t = R0*diag([J1,J2,J3])*R0';

% kiihtyvyysvoimavektorit FaccA FaccB
k1= (1-1/sqrt(3))/2;
k2= (1+1/sqrt(3))/2;
Yg1 = k2*Y1 + k1*Y2;
Yg2 = k1*Y1 + k2*Y2;
T1 =  MB_T(Yg1);
T2 =  MB_T(Yg2);

FaccB =  L/4*[Z3,I3,Z3,I3]'*T'*(R0JR0t*Tdot + mato(T*Ygdot)*R0JR0t*T)*[Z3,I3,Z3,I3]*[ddot(1:9); Y2dot];
M =  Arhoo*L/2*([k1*I3,Z3,k2*I3,Z3]'*[k1*I3,Z3,k2*I3,Z3] + [k2*I3,Z3,k1*I3,Z3]'*[k2*I3,Z3,k1*I3,Z3]) + L/2*([Z3,k1*I3,Z3,k2*I3]'*T2'*R0JR0t*T2*[Z3,k1*I3,Z3,k2*I3] + [Z3,k2*I3,Z3,k1*I3]'*T1'*R0JR0t*T1*[Z3,k2*I3,Z3,k1*I3]);
FaccA = M*[ddot2(1:9); Y2dot2];


if vaihto,
    % kinemaattinen matriisi, kun kiertymävektoria kutistetaan
    BB = [eye(9,9), zeros(9,3); zeros(3,9), Bvaihto];
    BBdot = [zeros(9,9), zeros(9,3); zeros(3,9), Bvaihtodot];
    FaccB = BB'*(FaccB + M*BBdot*ddot);
    M = BB'*M*BB;
    FaccA = M*ddot2;

    % disp('vaihto') % matlab-juttu
end    
