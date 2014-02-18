function Fii=BallConst(t,y0,v0,a0)

global MD;

  x0=[y0(1,1);y0(2,1);y0(3,1);y0(4,1);y0(5,1);y0(6,1)];
  Fii=x0;                               % tämän pitäisi olla oikein