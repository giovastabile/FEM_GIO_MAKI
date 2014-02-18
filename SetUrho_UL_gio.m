function Dummy=SetUrho_UL_gio(Urho)

global FEM;

Dummy=0;
% Asetetaan elementeille siirtymät
for Ele=1:FEM.NumEle
    FEM.ELE(Ele).UrhoE=zeros(12,1);
    for Node=1:FEM.ELE(Ele).NumNodes
        NodeJump = FEM.ELE(Ele).NDOF * (Node - 1);
        for dof=1:FEM.ELE(Ele).NDOF
            pos=(FEM.ID( FEM.ELE(Ele).Node(Node),dof ));
            if pos>0
                 FEM.ELE(Ele).UrhoE(NodeJump+dof)=Urho(pos);
             end
        end
    end
end