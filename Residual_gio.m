function r=Residual_gio(Urho)

Upda=0;

% residuaalivektorin r laskenta kts. kaava 2.85
Dummy=SetUrho_UL_gio(Urho);

% ulkoinen voimavektori
Fex = Fext_gio();

% sisäinen voimavektori
Fin = Fint_gio(Upda);

% residuaalivektori
r = Fex - Fin;
  
