function T = MB_T(y)
%
%	T = MB_T(y)
%
%	Laskee kiertym�liikkeen tangenttimatriisin T, kaava (4.11) 
%	
%	Input:
%		y 	kiertym�vektori, dim(y)=3
%
%   Output:
%      T    tangenttimatriisi, dim(T)=3x3
%

%
%	Tarkastettu: 5.11.2001, Jari M�kinen
%


UnitNormTreshold=1e-6;
y = y(:);

% kiertym�kulma
y_l = sqrt(y'*y);

ymato = [0,-y(3),y(2); y(3),0,-y(1); -y(2), y(1), 0];

% tangenttimatriisi T
if y_l < UnitNormTreshold;
    	% nollan ymp�rist�ss�
	T = eye(3) - 0.5*ymato + 1/6*ymato*ymato;
else
	T = sin(y_l)/y_l*eye(3,3) - (1 - cos(y_l))/y_l^2*ymato  + (y_l - sin(y_l))/y_l^3*y*y';
end
