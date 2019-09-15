function [referencia]=RefS_LMuslo(names, dinamica, LC_cadera, LC_rodilla)

%Hallar la posicion de la barra 1 y tomar los datos

X=Find_name(names,'l bar 1.X');
Y=Find_name(names,'l bar 1.Y');
Z=Find_name(names,'l bar 1.Z');

L_Bar1=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

eje_referencia={};

for x=1:length(L_Bar1(:,1))
    
    eje_ref = Eje_referencia(LC_rodilla(x,:),LC_cadera(x,:),L_Bar1(x,:));
    
    eje_referencia{x}=eje_ref;

      
end

referencia=eje_referencia;
  
end