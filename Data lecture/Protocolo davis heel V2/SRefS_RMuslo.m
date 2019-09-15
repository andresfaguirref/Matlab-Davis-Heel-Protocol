function [referencia]=SRefS_RMuslo(names, static, RC_cadera, RC_rodilla)

%Hallar la posicion de la barra 1 y tomar los datos

X=Find_name(names,'r bar 1.X');
Y=Find_name(names,'r bar 1.Y');
Z=Find_name(names,'r bar 1.Z');

R_Bar1=[static(:,X), static(:,Y), static(:,Z)];

eje_referencia={};

for x=1:length(R_Bar1(:,1))
    
    eje_ref = Eje_referencia(RC_rodilla(x,:),RC_cadera(x,:),R_Bar1(x,:));
    
    eje_referencia{x}=eje_ref;
    
end

referencia=eje_referencia;
  
end