function Dummy=SetUrho_UL(Urho,Ville,Adolf)

global FEM;

Dummy=0;  

% Place displacements inside elements

for Ele=1:FEM.NumEle
    FEM.ELE(Ele).UrhoE=zeros(12,1);
    for node=1:FEM.ELE(Ele).NumNodes
        NODE=FEM.ELE(Ele).Node(node);
        IDX_glob=FEM.ID(NODE,:);
        IDX_loc=(node-1)*6+1:1:6*node;
        FEM.ELE(Ele).UrhoE(IDX_loc)=Urho(IDX_glob);
        FEM.ELE(Ele).VilleE(IDX_loc)=Ville(IDX_glob);
        FEM.ELE(Ele).AdolfE(IDX_loc)=Adolf(IDX_glob);
    end
end
end

