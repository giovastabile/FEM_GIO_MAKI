function B=ConstraintG(t,y0,v0,a0)
% rajoiteyht�l�iden gradientti, pallorajoite heilurin alkup��ss�
global FEM;
  B=zeros(6,FEM.DOFS);
  B(1:6,1:6)=eye(6,6);