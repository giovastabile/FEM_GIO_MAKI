function MG=AddToGlobal(Ele,ML)
global FEM;
MG=zeros(FEM.DOFS);

for Node1=1:FEM.ELE(Ele).NumNodes
    Jump1=(Node1-1)*FEM.ELE(Ele).NDOF;
    for Dof1=1:FEM.ELE(Ele).NDOF
        for Node2=1:FEM.ELE(Ele).NumNodes
            for Dof2=1:FEM.ELE(Ele).NDOF
                Jump2=(Node2-1)*FEM.ELE(Ele).NDOF;
                IDR=FEM.ID ( FEM.ELE(Ele).Node(Node1),Dof1);
                IDC=FEM.ID( FEM.ELE(Ele).Node(Node2),Dof2);
                if IDR>0
                    if IDC>0
                        MG(IDR ,IDC)=ML((Node1-1)*Jump1+Dof1,(Node2-1)*Jump2+Dof2);
                    end
                end
            end
        end
    end
end

