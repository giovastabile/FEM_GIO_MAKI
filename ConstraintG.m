function B=ConstraintG(t,y0,v0,a0)
% rajoiteyhtälöiden gradientti, pallorajoite heilurin alkupäässä
global FEM;
  B=zeros(6,FEM.DOFS);
  B(1:6,1:6)=eye(6,6);