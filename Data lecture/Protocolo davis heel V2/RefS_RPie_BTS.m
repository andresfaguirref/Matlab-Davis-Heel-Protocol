function [referencia]=RefS_RPie_BTS(names, dinamic, RC_tobillo, RC_metatarso)

%Hallar la posicion de la barra 1 y tomar los datos

X=Find_name(names,'r heel.X');
Y=Find_name(names,'r heel.Y');
Z=Find_name(names,'r heel.Z');

R_Heel=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

X=Find_name(names,'r mall.X');
Y=Find_name(names,'r mall.Y');
Z=Find_name(names,'r mall.Z');

R_Mall=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

eje_referencia={};

for x=1:length(R_Heel(:,1))
    
    Vector_x=R_Heel(x,:)-RC_metatarso(x,:);
    Vector_z=RC_tobillo(x,:)-R_Mall(x,:);
    Vector_y=cross(Vector_z, Vector_x);
    
    Vector_z=cross(Vector_x, Vector_y);

    eje_referencia{x,1}=[Vector_x/norm(Vector_x); Vector_y/norm(Vector_y);Vector_z/norm(Vector_z)];
    
end

referencia=eje_referencia;
  
end