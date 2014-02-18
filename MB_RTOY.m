function Ys = MB_RTOY(Rs)
%   Ys = MB_logR(Rs)
%
%	Laskee kiertymävektorin kiertymämatriisista (rotaatiomatriisi), Rs = exp(mato(Ys))
%
% INPUT:
%   R    kiertymämatriisi, dim(R)=3x3
%
% OUTPUT:
%   Ys   kiertymävektori, dim(Ys) = 3  
%

% pieni kulma
fiiepsi = 1e-6;

% lasketaan kiertymämatriisia Rs vastaava kiertymävektori Ys, dim=3
% haetaan seuraavista muuttujista t1, t2, t3 suurin 
t1 = Rs(2,2)*Rs(3,3)-Rs(2,2)-Rs(3,3)+1-Rs(2,3)*Rs(3,2);
t2 = Rs(1,1)*Rs(3,3)-Rs(1,1)-Rs(3,3)+1-Rs(3,1)*Rs(1,3);
t3 = Rs(1,1)*Rs(2,2)-Rs(1,1)-Rs(2,2)+1-Rs(1,2)*Rs(2,1);
if t1 > t2,
    if t1 > t3,
        valinta = '1';
        tmax = t1;
    else
        valinta = '3';
        tmax = t3;
    end
else
    if t2 > t3,
        valinta = '2';
        tmax = t2;
    else
        valinta = '3';
        tmax = t3;
    end
end

if tmax < fiiepsi*fiiepsi,
    valinta = '4';
end  
  
switch  valinta
case '1',
    n(1) = tmax;
    n(2) = -Rs(1,2)*Rs(3,3)+Rs(1,2)+Rs(1,3)*Rs(3,2);
    n(3) = Rs(1,2)*Rs(2,3)-Rs(1,3)*Rs(2,2)+Rs(1,3);
    oppituus = 1/sqrt( n(1)*n(1) + n(2)*n(2) + n(3)*n(3) );
    n(1) = n(1)*oppituus;
    n(2) = n(2)*oppituus;
    n(3) = n(3)*oppituus;
case '2',    
    n(1) = -Rs(2,1)*Rs(3,3)+Rs(2,1)+Rs(2,3)*Rs(3,1);
    n(2) = tmax;
    n(3) = -Rs(1,1)*Rs(2,3)+Rs(2,3)+Rs(2,1)*Rs(1,3);
    oppituus = 1/sqrt( n(1)*n(1) + n(2)*n(2) + n(3)*n(3) );
    n(1) = n(1)*oppituus;
    n(2) = n(2)*oppituus;
    n(3) = n(3)*oppituus;
case '3',
    n(1) = Rs(3,2)*Rs(2,1)-Rs(2,2)*Rs(3,1)+Rs(3,1);
    n(2) = -Rs(1,1)*Rs(3,2)+Rs(3,2)+Rs(1,2)*Rs(3,1); 
    n(3) = tmax;
    oppituus = 1/sqrt( n(1)*n(1) + n(2)*n(2) + n(3)*n(3) );
    n(1) = n(1)*oppituus;
    n(2) = n(2)*oppituus;
    n(3) = n(3)*oppituus;
case '4',
    n(1) = 0;
    n(2) = 0;
    n(3) = 0;
end    
OP2 = 1/2;
cos_fiis = (Rs(1,1) + Rs(2,2) + Rs(3,3) -1) * OP2;
sin_fiis = (Rs(3,2)*n(1) - Rs(3,1)*n(2) - Rs(2,3)*n(1) + Rs(2,1)*n(3) + Rs(1,3)*n(2) - Rs(1,2)*n(3)) * OP2;
fiis = atan2(sin_fiis,cos_fiis);
Ys(1) = fiis*n(1);
Ys(2) = fiis*n(2);
Ys(3) = fiis*n(3);
% fiis voi olla myös negatiivinen, mikä otetaan huomioon seuraavasti 
%t1 = sign(fiis);  % merkkifunktio
%fiis = t1*fiis;
%sin_fiis = t1*sin_fiis;   % cosini on parillinen funktio

% pystyvektoroidaan
Ys = Ys(:);
