function Kt=Stiffness_UL(t,Urho,v0,a0,l0)

global FEM;

    Kt=zeros(FEM.DOFS,FEM.DOFS);
    Dummy=SetUrho_UL(Urho,v0,a0);
    for Ele=1:FEM.NumEle
        E=FEM.RC(FEM.ELE(Ele).Material).E;
        A=FEM.RC(FEM.ELE(Ele).Material).A;
        ElType=FEM.ELE(Ele).ElType;
        switch ElType
        case FEM.ET.ROD2D
          k_ele=Rod2D_K(FEM.ELE(Ele).X,FEM.ELE(Ele).UrhoE,E,A);
          MD=AddToGlobal(Ele,k_ele);
          Kt=Kt+MD;
        case FEM.ET.MAS3D    
        case FEM.ET.BEAM3D  
          G=.3*E;
          Mat(1)=E*A;
          Mat(2)=G*A;
          Mat(3)=G*A;
          Mat(4)=G*FEM.RC(FEM.ELE(Ele).Material).Iv;
          Mat(5)=E*FEM.RC(FEM.ELE(Ele).Material).Iy;
          Mat(6)=E*FEM.RC(FEM.ELE(Ele).Material).Iz;
          [Kmat,Kgeom]=B3DUL_Kint(Mat,FEM.ELE(Ele).X,FEM.ELE(Ele).S2,FEM.ELE(Ele).UrhoE,FEM.ELE(Ele).Rref,FEM.ELE(Ele).Kref);
          k_ele=Kmat+Kgeom;
          MD=AddToGlobal(Ele,k_ele);
          Kt=Kt+MD;
        end
    end
  
