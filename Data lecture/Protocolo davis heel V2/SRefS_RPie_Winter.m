function [referencia]=SRefS_RPie_Winter(names, static, RC_tobillo, RC_metatarso)

%Hallar la posicion de la barra 1 y tomar los datos

X=Find_name(names,'r heel.X');
Y=Find_name(names,'r heel.Y');
Z=Find_name(names,'r heel.Z');

R_Heel=[static(:,X), static(:,Y), static(:,Z)];

eje_referencia={};

for x=1:length(R_Heel(:,1))
    
    Vector_x=R_Heel(x,:)-RC_metatarso(x,:);
    Vector_int1=RC_tobillo(x,:)-R_Heel(x,:);
    Vector_int2=RC_metatarso(x,:)-R_Heel(x,:);
    
    Vector_z=cross(Vector_int1,Vector_int2);
    Vector_y=cross(Vector_z, Vector_x);

    eje_referencia{x,1}=[Vector_x/norm(Vector_x); Vector_y/norm(Vector_y);Vector_z/norm(Vector_z)];
    
end

referencia=eje_referencia;
  
end