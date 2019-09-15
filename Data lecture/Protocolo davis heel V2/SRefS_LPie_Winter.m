function [referencia]=SRefS_LPie_Winter(names, static, RC_tobillo, RC_metatarso)

%Hallar la posicion de la barra 1 y tomar los datos

X=Find_name(names,'l heel.X');
Y=Find_name(names,'l heel.Y');
Z=Find_name(names,'l heel.Z');

L_Heel=[static(:,X), static(:,Y), static(:,Z)];

eje_referencia={};

for x=1:length(L_Heel(:,1))
    
    Vector_x=L_Heel(x,:)-RC_metatarso(x,:);
    Vector_int1=RC_tobillo(x,:)-L_Heel(x,:);
    Vector_int2=RC_metatarso(x,:)-L_Heel(x,:);
    
    Vector_z=cross(Vector_int1,Vector_int2);
    Vector_y=cross(Vector_z, Vector_x);

    eje_referencia{x,1}=[Vector_x/norm(Vector_x); Vector_y/norm(Vector_y);Vector_z/norm(Vector_z)];
    
end

referencia=eje_referencia;
  
end