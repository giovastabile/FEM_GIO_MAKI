function [Ktot,Ctot,Kmat,Kgeom,Kcent,Cgyro] = B3D_KC(EA,GA2,GA3,GJ,EI2,EI3,Arhoo,J1,J2,J3,x0,suunta,d,ddot,ddot2)
%   [Ktot,Ctot,Kmat,Kgeom,Kcent,Cgyro] = B3D_KC(EA,GA2,GA3,GJ,EI2,EI3,Arhoo,J1,J2,J3,x0,suunta,d,ddot,ddot2)
%	 
%   Laskee Reissnerin 3D-palkin totaaliset jäykkyys- ja vaimennusmatriisit (Ktot, Ctot) 
%   käyttäen 1-pisteen gaussin integrointia, jossa on lineaariset interpolaatiofuktiot.  	
%	
%
% INPUT:
%   EA      veto-puristusjäykkyys, skalaari	
%   GA2     leikkausjäykkyys y-suuntaan, skalaari
%   GA3     leikkausjäykkyys z-suuntaan, skalaari
%   GJ      vääntöjäykkyys x-akselin ympäri, skalaari
%   EI2     taivutusjäykkyys y-akselin ympäri, skalaari
%   EI3     taivutusjäykkyys z-akselin ympäri, skalaari
%   Arhoo   A*rhoo = pituustiheys [kg/m], skalaari	
%   J1      hitausmomentti x-akselin ympäri, skalaari
%   J2      hitausmomentti y-akselin ympäri, skalaari
%   J3      hitausmomentti z-akselin ympäri, skalaari
%   x0      solmujen alkuasemavektori, dim(x0)=6, x0=[x01,y01,z01,x02,y02,z02]
%   suunta  paikallisen y-akselin suuntavektori, ei tarvitse olla ortogonaalinen palkin x-suunnan kanssa 
%   d       solmusiirtymävektori, dim(d)=12, d=[u1,v1,w1,fii_x1,fii_y1,fii_z1, u2,v2,w2,fii_x2,fii_y2,fii_z2] 
%   ddot    solmunopeusvektori, dim(ddot)=12
%   ddot2   solmukiihtyvyysvektori, dim(ddot2)=12
%
% OUTPUT:
%   Ktot    kokonaisjäykkyysmatriisi (sis. geom. mat. keskipakoismatriisit), dim(Ktot)=12x12
%   Ctot    kokonaisvaimennusmatriisi, dim(Ctot)=12x12
%   Kmat    materaalijäykkyysmatriisi, dim(Kmat)=12x12
%   Kgeom   geometrinen jäykkyysmatriisi, dim(Kgeom)=12x12
%   Kcent   keskipakoisjäykkyysmatriisi, dim(Kcent)=12x12
%   Cgyro   hyrrävaimennusmatriisi, dim(Cgyro)=12x12
%


%
%   tarkistettu 16.10.2002, Jari Mäkinen, <jari@mohr.me.tut.fi>
%

d = d(:); ddot = ddot(:); ddot2 = ddot2(:); suunta = suunta(:);
x0 = x0(:);
I3 = eye(3,3);
Z3 = zeros(3,3);

% pieni kulma: kun kulma on tätä pienempi lasketaan sisäinen voimavektori yksinkertaisemman kaavan mukaisesti
fii_epsi = 1e-6; 

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
        Bvaihto = (1-2*pi/fii_Y2)*I3+2*pi/fii_Y2^3*(Y2*Y2');
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
  

% palkin pituuden laskenta
L = norm(x0(4:6)-x0(1:3));
% palkin suuntainen yksikkövektori
E1 = 1/L*(-x0(1:3) + x0(4:6));

% kiertymäsuureet gaussin pisteessä
Yg = 1/2*(Y1 + Y2);
Ygdot = 1/2*(Y1dot + Y2dot);
Ygdot2 = 1/2*(Y1dot2 + Y2dot2);
% kiertymäsuureiden derivaatat gaussin pisteessä
Yg_s = 1/L*(Y2 - Y1); 

% rotaatiomatriisi R, transformaatiomatriisit T 
R = MB_R(Yg);
T = MB_T(Yg);
Tdot = MB_Tdot(Ygdot,Yg);

% kimmoviivan derivaatta gaussin pisteessä 
x0dp_dot = 1/L*(x0(4:6)+d(7:9)-x0(1:3)-d(1:3)); 
% vektori transpoosi(R)*xdot
Rtxdot = R'*x0dp_dot;
% venymävektori GAMMA ;
G = Rtxdot - E1;
% Kaareutumavektori K (matriisi T kertaa vektori Yg_s
K = T*Yg_s;

% kimmomatriisi R0ER0t (R^6x6) ja hitausmatriisi (R^3x3) 
E2 = suunta - E1'*suunta*E1;
E2 = E2/norm(E2);
E3 = cross(E1,E2);
R0 = [E1,E2,E3];
R0JR0t = R0*diag([J1,J2,J3])*R0';
R0ER0t = [R0*diag([EA,GA2,GA3])*R0', zeros(3,3); zeros(3,3), R0*diag([GJ,EI2,EI3])*R0'];

% sisäinen normaalivoima ja momentti NjaM ja jäykkyysmatriisi
N = R0*diag([EA,GA2,GA3])*R0'*G;
M = R0*diag([GJ,EI2,EI3])*R0'*K;
C1 = MB_C1(Yg_s,Yg);
Bh = [R', Z3, mato(Rtxdot)*T; Z3, T, C1];
Lambda = [-1/L*eye(6,6), 1/L*eye(6,6); Z3, 1/2*I3, Z3, 1/2*I3];
B = Bh*Lambda;
C2_1 = MB_C2(M,Yg);
C2_2 = MB_C2(mato(N)*Rtxdot,Yg);
C3 = MB_C3(M,Yg_s,Yg); 
% materiaalinen jäykkyysmatriisi
Kmat = L*B'*R0ER0t*B;
% geometrinen jäykkyyematriisi
Kg = [Z3,Z3,-R*mato(N)*T; Z3,Z3,C2_1; T'*mato(N)*R', C2_1', C3+C2_2+T'*mato(N)*mato(Rtxdot)*T]; 
Kgeom = L*Lambda'*Kg*Lambda;

% hitausvoimavektorit ja tangenttimatriisit
Omega = T*Ygdot;
A = Tdot*Ygdot;
C2 = MB_C2(mato(Omega)*R0JR0t*Omega+R0JR0t*A, Yg);
% kiihtyvyysvoimavektorit FaccA FaccB
k1= (1-1/sqrt(3))/2;
k2= (1+1/sqrt(3))/2;
Yg1 = k2*Y1 + k1*Y2;
Yg2 = k1*Y1 + k2*Y2;
Yg1dot2 = k2*Y1dot2 + k1*Y2dot2;
Yg2dot2 = k1*Y1dot2 + k2*Y2dot2;
T1 =  MB_T(Yg1);
T2 =  MB_T(Yg2);

Kcent = L/4*[Z3,I3,Z3,I3]'*(C2+T'*((mato(Omega)*R0JR0t-mato(R0JR0t*Omega))*MB_C1(Ygdot,Yg) + R0JR0t*(MB_C5(Ygdot,Yg)+0*MB_C1(Ygdot2,Yg))))*[Z3,I3,Z3,I3] +   L/2*([Z3,k1*I3,Z3,k2*I3]'*(MB_C2(R0JR0t*T2*Yg2dot2, Yg2) + T2'*R0JR0t*MB_C1(Yg2dot2,Yg2))*[Z3,k1*I3,Z3,k2*I3] + [Z3,k2*I3,Z3,k1*I3]'*(MB_C2(R0JR0t*T1*Yg1dot2, Yg1) + T1'*R0JR0t*MB_C1(Yg1dot2,Yg1))*[Z3,k2*I3,Z3,k1*I3]); 
Cgyro = L/4*[Z3,I3,Z3,I3]'*T'*((mato(Omega)*R0JR0t-mato(R0JR0t*Omega))*T + R0JR0t*(Tdot+MB_C1(Ygdot,Yg)))*[Z3,I3,Z3,I3];

if vaihto,
    % disp('vaihto')
    % sisäinen voimavektori 
    [FaccA, FaccB, M] = B3D_Facc(Arhoo,J1,J2,J3,x0,suunta,[d(1:9); Y2],[ddot(1:9); Y2dot],[ddot2(1:9); Y2dot2]); 
    % hitausvoimavektori
  
    F_int = B3D_Fint(EA,GA2,GA3,GJ,EI2,EI3,x0,suunta,[d(1:9); Y2]);
    Facc = FaccA + FaccB;

    Y2 = d(10:12); 
    Y2dot = ddot(10:12);
    Y2dot2 = ddot2(10:12);

    fii2 = norm(Y2);
    e = Y2/fii2;
    fii2dot = 1/fii2*(Y2dot'*Y2); 
    edot = 1/fii2*Y2dot - 1/fii2^3*(Y2'*Y2dot)*Y2;

    % kinemaattinen matriisi, kun kiertymävektoria kutistetaan
    %BB = [ I3, Z3; Z3, Bvaihto];
    Z99 = zeros(9,9); Z93 = zeros(9,3);
    BB = [eye(9,9), Z93; Z93', Bvaihto];
    BBdot = [Z99, Z93; Z93', Bvaihtodot];
    BBdot2 = [Z99, Z93; Z93', 2*pi/fii2^2*( Y2dot2*e' + Y2dot*edot' + e*Y2dot2' + edot*Y2dot' - 3*(e'*Y2dot2 + edot'*Y2dot)*e*e' - 3*(e'*Y2dot)*(edot*e' + e*edot') + (e'*Y2dot2 + edot'*Y2dot)*I3) - 2*fii2dot/fii2*Bvaihtodot];
      
    dBF = 2*pi/fii2^2*( Facc(10:12)*e' + e*Facc(10:12)' - 3*(Facc(10:12)'*e)*(e*e') + (e'*Facc(10:12))*I3 ); 
    Kcent = BB'*( Kcent*BB + Cgyro*BBdot + M*BBdot2) + [Z99, Z93; Z93', dBF];
    Cgyro = BB'*(Cgyro*BB + 2*M*BBdot);
        %%% voidaan myös yhdistää: 
        %%%  F = FaccA + FaccB + F_int;
        %%% Kcent = BB'*( Kcent*BB + Cgyro*BBdot + M*BBdot2); 
        %%% Kg1= 2*pi/fii_Y2^2*( F(10:12)*e' + e*F(10:12)' - 3*(F(10:12)'*e)*(e*e') + (e'*F(10:12))*I3);
    
    % geometrinen jäykkyysmatriisi
    Kg = BB'*Kgeom*BB;
    e = Y2/norm(Y2);
    Kg1= 2*pi/fii_Y2^2*( F_int(10:12)*e' + e*F_int(10:12)' - 3*(F_int(10:12)'*e)*(e*e') + (e'*F_int(10:12))*I3);
    Kgeom = Kg + [Z99, Z93; Z93', Kg1];
    %materiaalinen jäykkyysmatriisi    
    Kmat = BB'*Kmat*BB;
end    

Ktot = Kmat + Kgeom + Kcent;
Ctot = Cgyro;

