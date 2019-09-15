function longitud=LStep_length(names, dinamic, hs1, hs2)


X=Find_name(names,'l heel.X');
Y=Find_name(names,'l heel.Y');
Z=Find_name(names,'l heel.Z');

L_Heel=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

primero=L_Heel(hs1,:);
segundo=L_Heel(hs2,:);

vector=segundo-primero;
longitud=norm(vector);
end