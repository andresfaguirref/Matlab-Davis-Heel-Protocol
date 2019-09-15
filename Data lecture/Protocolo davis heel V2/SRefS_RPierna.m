function [referencia]=SRefS_RPierna(names, static, RC_rodilla, RC_tobillo)

%Hallar la posicion de la barra 1 y tomar los datos

X=Find_name(names,'r bar 2.X');
Y=Find_name(names,'r bar 2.Y');
Z=Find_name(names,'r bar 2.Z');

R_Bar2=[static(:,X), static(:,Y), static(:,Z)];

eje_referencia={};

for x=1:length(R_Bar2(:,1))
    
    Vector_x=RC_rodilla(x,:)-RC_tobillo(x,:);
    Vector_int1=R_Bar2(x,:)-RC_rodilla(x,:);
    Vector_int2=RC_tobillo(x,:)-RC_rodilla(x,:);
    
    
    Vector_y=cross(Vector_int1,Vector_int2);
    Vector_z=cross(Vector_x, Vector_y);

    eje_referencia{x,1}=[Vector_x/norm(Vector_x); Vector_y/norm(Vector_y);Vector_z/norm(Vector_z)];
    
end

referencia=eje_referencia;
  
end