function r=Residual_UL_prova(Urho)
t=0;
Ville=Urho;
Adolf=Urho;
global FEM
Upda=0;
Lam=[0;0;0];

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
r = FEM.EXloads+Fex - Fin - B'*Lam;%
  
