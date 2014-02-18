function Tdot = Tdot(ydot,y)
%
%	Tdot=MB_Tdot(ydot,y)
%	Laskee tangenttioperaattorin  T aikaderivaatan
%
%	Input:
%		ydot	kiertymävektorin aikaderivaatta		
%		y 	kiertymävektori
%
%	Output:
%		Tdot	T:n aikaderivaatta
%

%
%	Tarkastettu: 5.11.2001, Jari Mäkinen
%


ydot=ydot(:);
y=y(:);

w=sqrt(y'*y);

if w < 1e-6,
  Tdot = -1/2*MB_mato(ydot);
else
  cw = cos(w);
  sw = sin(w);
  c1 = (w*cw-sw)/w^3;
  c2 = (w*sw+2*cw-2)/w^4;
  c3 = (3*sw-2*w-w*cw)/w^5;
  c4 = (cw-1)/w^2;
  c5 = (w-sw)/w^3;
  ydoty = ydot'*y;
  Tdot = c1*ydoty*eye(3,3) - c2*ydoty*MB_mato(y) + c3*ydoty*y*y' + c4*MB_mato(ydot) + c5*(ydot*y' + y*ydot');
end
 
