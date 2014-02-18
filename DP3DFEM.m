% DP3DFEM

%Heilurin ratkaisu kolmiulotteisena
%k‰ytt‰en palkkielementti‰ ja pistemassaa solmussa 2. Hessu 21.4.05

% Globaalimuuttuja FEM siirt‰‰ mallin geometria ja topologiatiedon aliohjelmille
% clear all;
clear FEM;
global FEM;

% elementtityypit
%
% FEM alustus
FEM.PARA.g=9.82;
FEM.PARA.GRAVITY=1;
FEM.PARA.Gx=0;
FEM.PARA.Gy=-9.82;
FEM.PARA.Gz=0;
FEM.ET.ROD2D=1;
FEM.ET.MAS3D=2;
FEM.ET.BEAM3D=3;

% Section properties
B=.008;
H=.008;
Ls=0.4; %Pendulum Length

FEM.RC(1).E=2.E11;     % Young's Modulus
FEM.RC(1).A=B*H;       % Area of the cross section 
FEM.RC(1).Iz=B*H^3/12; % Moment of inertia respect to  
FEM.RC(1).Iy=B^3*H/12;
FEM.RC(1).Iv=.14*B^2*H;
FEM.RC(1).rho=.1/Ls; % sauvan pituusmassa kg/m
FEM.RC(1).m=1; % pistemassa kg
% Hitausmatriisi
Roo=.1/FEM.RC(1).A/Ls;                  % pakotettu siten, ett‰ massa on 0.1 kg
FEM.RC(1).J=GetBeamJ(Roo,FEM.RC(1).Iy,FEM.RC(1).Iz);
% solmukoordinaatit
% Elementit
FEM.NumEle=2;                                       % Sauvamainen palkki ja pistemassa
FEM.NumNodes=2;

% Model parameters
% ID-Tables
for i=1:6
  FEM.ID(1,i)=i;
  FEM.ID(2,i)=i+6;
end

FEM.DOFS=12;                                         % mallin vapausasteiden lukum‰‰r‰

% Toisen solmun alkuasema
Alfa=0;                                             % X-akselin ymp‰ri
Beta=0;                                             % Y-akselin ymp‰ri
Gamma=-2*pi/180;                                   % Z-akselin ymp‰ri

% GetBaseVectors

Y=[Alfa,0,0];
Rx=MB_R(Y);

Y=[0,Beta,0];
Ry=MB_R(Y);

Y=[0,0,Gamma];
Rz=MB_R(Y);

R=Rz*Ry*Rx;

E=[0,1,0]';
FEM.ELE(1).S2=R*E;

% Alkutilan solmukahvelit

for i=1:FEM.NumNodes
  FEM.Kahveli(i).Y=[0,0,0]';
end

% Solmukoordinaatit
Node(1).x=0;
Node(1).y=0;
Node(1).z=0;

temp=R*[Ls,0,0]';
Node(2).x=temp(1,1);
Node(2).y=temp(2,1);
Node(2).z=temp(3,1);


% Palkkielementti
  FEM.ELE(1).NumNodes=2;
  FEM.ELE(1).NDOF=6;
  FEM.ELE(1).Node(1)=1;
  FEM.ELE(1).Node(2)=2;
  FEM.ELE(1).Material=1;
  FEM.ELE(1).ElType=3;      % elementtityyppi palkki 3D
  FEM.ELE(1).Rref=eye(3);
  FEM.ELE(1).Kref=zeros(3,1);

  for i=1:1
    x1=Node(FEM.ELE(i).Node(1)).x;
    y1=Node(FEM.ELE(i).Node(1)).y;
    z1=Node(FEM.ELE(i).Node(1)).z;
    x2=Node(FEM.ELE(i).Node(2)).x;
    y2=Node(FEM.ELE(i).Node(2)).y;
    z2=Node(FEM.ELE(i).Node(2)).z;
    FEM.ELE(i).X=[x1,y1,z1,x2,y2,z2];
  end
  
% Pistemassa
  FEM.ELE(2).NumNodes=1;
  FEM.ELE(2).NDOF=6;
  FEM.ELE(2).Node(1)=2;
  FEM.ELE(2).Material=1;
  FEM.ELE(2).ElType=2;      % elementtityyppi palkki 3D
  FEM.ELE(2).Rref=eye(3);

  for i=1:FEM.NumNodes
    FEM.Kahveli(i).R=eye(3);  
  end
  
% Aikaintegrointiarvot

t0=0;
t_max=1.5;
h=.005;                     % aika-askel
tspan=[t0;t_max;h];         % vector

% Alkuehdot

[y0,v0,l0]=GetInitialState;

[tout,yout,vout,rout,lout,stats] = NewmarkC3D('Residual_UL','TangentMatrices','BallConst',tspan,y0,v0,l0);
plot(tout,yout(:,6))

