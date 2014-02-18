function M=Mass(t,y0,v0,a0)

global FEM;
% Massamatriisi
    M=zeros(FEM.DOFS);
    for Ele=1:FEM.NumEle
      EleType=FEM.ELE(Ele).ElType;
      switch EleType
        case FEM.ET.ROD2D
          m_ele=Rod2D_M(FEM.ELE(Ele).X,FEM.RC( FEM.ELE(Ele).Material).rho);
        case FEM.ET.MAS3D  
          m_ele=Mass3D_M(FEM.RC( FEM.ELE(Ele).Material).m);
        case FEM.ET.BEAM3D 
          Mat(1)=FEM.RC( FEM.ELE(Ele).Material).rho ;  
          Mat(2)=FEM.RC( FEM.ELE(Ele).Material).J(1,1);  
          Mat(3)=FEM.RC( FEM.ELE(Ele).Material).J(2,2);  
          Mat(4)=FEM.RC( FEM.ELE(Ele).Material).J(3,3);  
          m_ele=B3DUL_M(Mat,FEM.ELE(Ele).X,FEM.ELE(Ele).S2,FEM.ELE(Ele).UrhoE,FEM.ELE(Ele).VilleE,FEM.ELE(Ele).AdolfE);
      end % switch

      MG=AddToGlobal(Ele,m_ele);
      M=M+MG;
    end

    % end function