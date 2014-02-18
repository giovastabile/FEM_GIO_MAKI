function Fin=Finer(t,y0,v0,a0)

global FEM;
% sisäisen voimavektorin laskenta
Fin=zeros(FEM.DOFS,1);

    for Ele=1:FEM.NumEle

        E=FEM.RC(FEM.ELE(Ele).Material).E;
        A=FEM.RC(FEM.ELE(Ele).Material).A;

        switch FEM.ELE(Ele).ElType
        case FEM.ET.MAS3D
          f_ele=Mass3D_M(FEM.RC( FEM.ELE(Ele).Material).m)*a0(7:12);  
        case FEM.ET.BEAM3D
          Mat(1)=FEM.RC( FEM.ELE(Ele).Material).rho ;  
          Mat(2)=FEM.RC( FEM.ELE(Ele).Material).J(1,1);  
          Mat(3)=FEM.RC( FEM.ELE(Ele).Material).J(2,2);  
          Mat(4)=FEM.RC( FEM.ELE(Ele).Material).J(3,3);  
          f_ele=B3DUL_Finer(Mat,FEM.ELE(Ele).X,FEM.ELE(Ele).S2,FEM.ELE(Ele).UrhoE,FEM.ELE(Ele).VilleE,FEM.ELE(Ele).AdolfE);
           
         end % switch    

         Jump=FEM.ELE(Ele).NDOF;
          for Node=1:FEM.ELE(Ele).NumNodes
            for Dof1=1:FEM.ELE(Ele).NDOF
                IDX=FEM.ID(FEM.ELE(Ele).Node(Node),Dof1);
                if IDX>0
                  Fin(IDX)=Fin(IDX)+f_ele((Node-1)*Jump+Dof1);
                end
            end % for
          end % for
    
          
            
    end % for
    
    
    % end function

