function Fint=Fint_gio(UpdateReference)

global FEM;
% sisäisen voimavektorin laskenta
Fint=zeros(FEM.DOFS,1);

    for Ele=1:FEM.NumEle

        E=FEM.RC(FEM.ELE(Ele).Material).E;
        A=FEM.RC(FEM.ELE(Ele).Material).A;

        switch FEM.ELE(Ele).ElType
        case FEM.ET.ROD2D
          f_ele=Rod2D_F(FEM.ELE(Ele).X,FEM.ELE(Ele).UrhoE,E,A);
        case 4 %FEM.ET.Beam3D
          G=.3*E;
          EIz=E*FEM.RC(FEM.ELE(Ele).Material).Iz;
          EIy=E*FEM.RC(FEM.ELE(Ele).Material).Iy;
          Iv=FEM.RC(FEM.ELE(Ele).Material).Iv;
          f_ele=B3D_Fint(G*A,G*A,G*A,G*Iv,EIz,EIy,FEM.ELE(Ele).X,S2,FEM.ELE(Ele).UrhoE);
        case FEM.ET.MAS3D
          f_ele=zeros(6,1);  
        case FEM.ET.BEAM3D
          G=.3*E;
          Mat(1)=E*A;
          Mat(2)=G*A;
          Mat(3)=G*A;
          Iv=FEM.RC(FEM.ELE(Ele).Material).Iv;
          Mat(4)=G*Iv;
          Mat(5)=E*FEM.RC(FEM.ELE(Ele).Material).Iy;
          Mat(6)=E*FEM.RC(FEM.ELE(Ele).Material).Iz;
          [f_ele,Rref,Kref]=B3DUL_Fint(Mat,FEM.ELE(Ele).X,FEM.ELE(Ele).S2,FEM.ELE(Ele).UrhoE,FEM.ELE(Ele).Rref,FEM.ELE(Ele).Kref);

          if UpdateReference>0
              FEM.ELE(Ele).Rref=Rref;
              FEM.ELE(Ele).Kref=Kref;
          end    
            
        end    
            
        Jump=FEM.ELE(Ele).NDOF;
        
        for node=1:FEM.ELE(Ele).NumNodes
            NODE=FEM.ELE(Ele).Node(node);
            IDX_glob=FEM.ID(NODE,:);
            IDX_loc=(node-1)*6+1:1:6*node;
            Fint(IDX_glob)=Fint(IDX_glob)+f_ele(IDX_loc);
        end
            
        
        
%         for Node=1:FEM.ELE(Ele).NumNodes
%           for Dof1=1:FEM.ELE(Ele).NDOF
%               IDX=FEM.ID(FEM.ELE(Ele).Node(Node),Dof1);
%               if IDX>0
%               Fint(IDX)=Fint(IDX)+f_ele((Node-1)*Jump+Dof1);
%               end
%           end
%         end
    end
    
    
    

