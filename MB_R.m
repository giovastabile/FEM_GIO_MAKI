function R = MB_R(y)
%
%	R = MB_R(y)
%
%	Laskee kiertymismatriisin R, miss� y on kiertym�vektori, kaava (4.5)
%
%	Input:
%		y 	kiertym�vektori, dim(y)=3
%
%   Output:
%      R    kiertymismatriisi, dim(R)=3x3
%

%
%	Tarkastettu: 5.11.2001, Jari M�kinen
%

UnitNormTreshold=1e-6;
y=y(:);

% vinosymmetrinen matriisi (crossprod(y,a) = ymato*a)
ymato=[0,-y(3),y(2); y(3),0,-y(1); -y(2), y(1), 0];

% kiertymiskulma
y_l = sqrt(y'*y);

% kiertymismatriisi 
if y_l < UnitNormTreshold,
    % nollan ymp�rist�ss�
	R = eye(3,3) + ymato + 0.5*ymato*ymato;
else
	R = eye(3,3) + sin(y_l)/y_l*ymato + (1-cos(y_l))/y_l^2*ymato*ymato;
end