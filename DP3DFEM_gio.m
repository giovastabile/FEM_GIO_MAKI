% Geometrically exact Simo-Reissner Beam

clear all;
clear FEM;
global FEM;

% Element type and properties

%% FEM parameters

FEM.PARA.g=0;           % Gravity
FEM.PARA.GRAVITY=1;     % Gravity factor
FEM.PARA.Gx=0;          % Gravity in x-Direction
FEM.PARA.Gy=-0;         % Gravity in y-Direction
FEM.PARA.Gz=0;          % Gravity in y-Direction
FEM.ET.ROD2D=1;         % Rod 2D
FEM.ET.MAS3D=2;         % Mass3D
FEM.ET.BEAM3D=3;        % Beam3D

%% Section properties

B=.05;                 % Height Section
H=.05;                 % Width Section
Ls=10;                 % Lenght of the Beam

FEM.RC(1).E=2.E11;      % Young's Modulus
FEM.RC(1).A=B*H;        % Area of the cross section
FEM.RC(1).Iz=B*H^3/12;  % Moment of inertia respect to z
FEM.RC(1).Iy=B^3*H/12;  % Moment of inertia respect to y
FEM.RC(1).Iv=.14*B^2*H; % Torsional moment of Inertia
FEM.RC(1).rho=.1/Ls;    % Mass of the Beam per unit of length kg/m
FEM.RC(1).m=1;          % Concentrated Mass kg

%% Inertia beam properties

Roo=.1/FEM.RC(1).A/Ls;                                  % Scalar Mass Density
FEM.RC(1).J=GetBeamJ(Roo,FEM.RC(1).Iy,FEM.RC(1).Iz);    % Inertia tensor
%% Fem nodes coordinates

% Number of Elements and Number of Nodes
FEM.NumEle=5;                                       
FEM.NumNodes=6;

%% Model parameters

% ID-Tables, It gives the global ID of each dof and each node.

for i=1:FEM.NumNodes;
    FEM.ID(i,:)=(6*(i-1)+1:6*i);
end

FEM.DOFS=FEM.NumNodes*6;        % DOF number
%% To use only of the pendulum example.
%Remove the block comment to use.
%{

% Starting Position of the second Node of the beam.
Alfa=0;                         % Angle respect to X
Beta=0;                         % Angle respect to Y
Gamma=-2*pi/180;                % Angle respect to Z

% GetBaseVectors

Y=[Alfa,0,0];
Rx=MB_R(Y);

Y=[0,Beta,0];
Ry=MB_R(Y);

Y=[0,0,Gamma];
Rz=MB_R(Y);

R=Rz*Ry*Rx;

E=[0,1,0]';
for i=1:FEM.NumEle
    FEM.ELE(i).S2=R*E;
end

% Initial Condition

for i=1:FEM.NumNodes
    FEM.Kahveli(i).Y=[0,0,0]';
end

%}
%% Node coordinates

for i=1:FEM.NumNodes
    
    Node(i).x=0+(i-1)*Ls/(FEM.NumNodes-1);
    Node(i).y=0;
    Node(i).z=0;
end

% Node(1).x=0;
% Node(1).y=0;
% Node(1).z=0;
% 
% Node(2).x=Ls/3;
% Node(2).y=0;
% Node(2).z=0;
% 
% Node(3).x=2*Ls/3;
% Node(3).y=0;
% Node(3).z=0;
% 
% Node(4).x=Ls;
% Node(4).y=0;
% Node(4).z=0;

%% Beam Elements

for i=1:FEM.NumEle
    FEM.ELE(i).NumNodes=2;
    FEM.ELE(i).NDOF=6;
    FEM.ELE(i).Node(1)=i;
    FEM.ELE(i).Node(2)=i+1;
    FEM.ELE(i).Material=1;
    FEM.ELE(i).ElType=3;            % Element Type
    FEM.ELE(i).Rref=eye(3);
    FEM.ELE(i).Kref=zeros(3,1);
    FEM.ELE(i).S2=[0 1 0];
end


for i=1:FEM.NumEle;
    x1=Node(FEM.ELE(i).Node(1)).x;
    y1=Node(FEM.ELE(i).Node(1)).y;
    z1=Node(FEM.ELE(i).Node(1)).z;
    x2=Node(FEM.ELE(i).Node(2)).x;
    y2=Node(FEM.ELE(i).Node(2)).y;
    z2=Node(FEM.ELE(i).Node(2)).z;
    FEM.ELE(i).X=[x1,y1,z1,x2,y2,z2];
end

%% Concentrated mass elements

% Pistemassa
%   FEM.ELE(2).NumNodes=1;
%   FEM.ELE(2).NDOF=6;
%   FEM.ELE(2).Node(1)=2;
%   FEM.ELE(2).Material=1;
%   FEM.ELE(2).ElType=2;      % Element Type
%   FEM.ELE(2).Rref=eye(3);
%% Set external loads

% Insert external loads FEM.Loads(node,:)=[N Ty Tz Mt My Mz] 

FEM.Loads=zeros(FEM.NumNodes,6);
FEM.Loads(6,:)=[0 0 0 0 0 1];
FEM.EXloads=zeros(36,80001);



for t=1:80001;
for node=1:FEM.NumNodes
    ID_glob=FEM.ID(node,:)';
    FEM.EXloads(ID_glob,t)=65450/80001*t*FEM.Loads(node,:);
end
end




%% Set external Costraints

% Insert costraints for each node for each degree of freedom.
% Fem.Costraints(node,:)=[TRx TRy TRz Rx Ry Rz] 1 for free 0 for blocked

FEM.Costraints=ones(FEM.NumNodes,6);
FEM.Costraints(1,:)=[0 0 0 0 0 0]; % Node 1 is blocked




%%

for i=1:FEM.NumNodes
    FEM.Kahveli(i).R=eye(3);
end

% Time Parameters

t0=0;                       % Initial time
t_max=400;                  % Final time
h=.005;                     % Time increment
tspan=[t0;t_max;h];         % Paramaters for time vector


% Initial conditions

[y0,v0,l0]=GetInitialState(0);


% options=optimset('Display','iter','MaxFunEvals',200000,'MaxIter',2000,'TolX',1e-6);
% eta=fsolve(@Residual_UL_prova,y0,options);

%% Solve Using Newmark method

[tout,yout,vout,rout,lout,stats] = NewmarkC3D('Residual_UL','TangentMatrices','BallConst',tspan,y0,v0,l0);

% plot(tout,yout(:,62))

