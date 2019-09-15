function longitud=RStep_length(names, dinamic, hs1, hs2)


X=Find_name(names,'r heel.X');
Y=Find_name(names,'r heel.Y');
Z=Find_name(names,'r heel.Z');

R_Heel=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

primero=R_Heel(hs1,:);
segundo=R_Heel(hs2,:);

vector=segundo-primero;
longitud=norm(vector);
end