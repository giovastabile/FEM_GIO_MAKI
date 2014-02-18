function r=Residual(t,Urho,Ville,Adolf,Lam,Upda)

global FEM
global Counter

% Set the displacement and the acceleration of each element
Dummy=SetUrho_UL(Urho,Ville,Adolf);

% External Force Vector
Fex = Fext(t,Urho);

% Internal force Vector
Fin = Fint(t,Urho,Ville,Adolf,Upda);

% hitausvoimavektori
%Finert = Finer(t,Urho,Ville,Adolf);

% sidosehtojen kinemaattinen matriisi
B = ConstraintG(t,Urho,Ville,Adolf);

% residuaalivektori
r = FEM.EXloads(:,Counter)+Fex - Fin - B'*Lam;%-Finert;
  
