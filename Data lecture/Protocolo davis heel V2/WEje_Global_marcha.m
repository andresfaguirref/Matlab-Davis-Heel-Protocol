function eje=WEje_Global_marcha(names,dinamic, rhs1, rhs2, lhs1, lhs2)

Eje_Global=[1 0 0; 0 1 0; 0 0 1];

X=Find_name(names,'sacrum.X');
Y=Find_name(names,'sacrum.Y');
Z=Find_name(names,'sacrum.Z');

Sacrum=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

Evento_RHS=[Sacrum(rhs1,:);Sacrum(rhs2,:)];
Evento_LHS=[Sacrum(lhs1,:);Sacrum(lhs2,:)];

vector_r=Evento_RHS(2,:)-Evento_RHS(1,:);
vector_r=vector_r/norm(vector_r);

vector_l=Evento_LHS(2,:)-Evento_LHS(1,:);
vector_l=vector_l/norm(vector_l);


vector=vector_r+vector_l;
vector=vector/norm(vector);

vGAITK=cross(Eje_Global(2,:),vector);
vGAITJ=cross(vGAITK,Eje_Global(2,:));

eje=[Eje_Global(2,:); vGAITJ; vGAITK];


end