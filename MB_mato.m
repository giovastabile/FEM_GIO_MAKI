function Y=MB_mato(y)
%
%	Y = MB_mato(y)
%	Muuttaa vektorin y vinosymmetriseksi matriisiksi Y, 
%	s.e. Ya=y x a, missä a on mielivaltainen vektori
%
%
%
%
%

Y=[0,-y(3),y(2); y(3),0,-y(1); -y(2), y(1), 0];