function C=Damping(t,y0,v0,a0)

global FEM;
  C=zeros(FEM.DOFS,FEM.DOFS);