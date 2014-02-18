function Fex=Fext(t,Urho)

global FEM;

Fex=zeros(FEM.DOFS,1);

    if FEM.PARA.GRAVITY==1
    for Ele=1:FEM.NumEle
        rho=FEM.RC(FEM.ELE(Ele).Material).rho;
        switch FEM.ELE(Ele).ElType
          case FEM.ET.BEAM3D  
            FEM.ELE(1,Ele).fg=Beam3D_G(FEM.ELE(Ele).X,rho);
          case FEM.ET.MAS3D
            FEM.ELE(1,Ele).fg(1,1)=FEM.RC(FEM.ELE(Ele).Material).m*FEM.PARA.Gx;  
            FEM.ELE(1,Ele).fg(2,1)=FEM.RC(FEM.ELE(Ele).Material).m*FEM.PARA.Gy;  
            FEM.ELE(1,Ele).fg(3,1)=FEM.RC(FEM.ELE(Ele).Material).m*FEM.PARA.Gz;  
            FEM.ELE(1,Ele).fg(4,1)=0;  
            FEM.ELE(1,Ele).fg(5,1)=0;  
            FEM.ELE(1,Ele).fg(6,1)=0;  
        end %switch  
        Jump=FEM.ELE(Ele).NDOF;
        for Node=1:FEM.ELE(Ele).NumNodes
          for Dof1=1:FEM.ELE(Ele).NDOF
              IDX=FEM.ID(FEM.ELE(Ele).Node(Node),Dof1);
              if IDX>0
              FEM.ELE(1,Ele).Fex(IDX)=Fex(IDX)+FEM.ELE(1,Ele).fg((Node-1)*Jump+Dof1);
              end
          end
        end
    end
    end
for i=1:FEM.NumEle
    FEM.ELE(1,i).fg=[zeros(6*(i-1),1);FEM.ELE(1,i).fg;zeros((FEM.DOFS-6*(i-1)-12),1)];
end

for i=1:FEM.NumEle
    Fex=Fex+FEM.ELE(1,i).fg;
end
    
    
   