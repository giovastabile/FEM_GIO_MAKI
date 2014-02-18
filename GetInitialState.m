%   GetInitialState
function [q0,qdot0,l0]=GetInitialState(par)

global FEM
a=FEM.Costraints';
if nargin<1
    par=0.0001;
end
% initial position vector and initial velocity vector

q0=a(:)*par;
qdot0=a(:)*par;

% Lagrange multipliers

for i=1:6
    l0(i,1)=0;
end

% end function