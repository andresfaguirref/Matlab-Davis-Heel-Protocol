function [referencia, mid_asis]=RefS_Pelvis(names, dinamica)

%Hallar el sistema de referencia del tronco

X=Find_name(names,'r asis.X');
Y=Find_name(names,'r asis.Y');
Z=Find_name(names,'r asis.Z');

R_Asis=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l asis.X');
Y=Find_name(names,'l asis.Y');
Z=Find_name(names,'l asis.Z');

L_Asis=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'sacrum.X');
Y=Find_name(names,'sacrum.Y');
Z=Find_name(names,'sacrum.Z');

Sacrum=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

punto_nuevo=[];
vector=L_Asis-R_Asis;


for x=1:length(L_Asis(:,1))  
 
    distancia=norm(vector(x,:));
   
    punto_nuevo(x,:)=Mover_punto(vector(x,:),R_Asis(x,:), distancia/2); 
    
    
end

mid_asis=punto_nuevo;

eje_referencia={};

for x=1:length(L_Asis(:,1))
    
    Vector_z=L_Asis(x,:)-punto_nuevo(x,:);
    Vector_int=Sacrum(x,:)-punto_nuevo(x,:);

    Vector_x=cross(Vector_z,Vector_int);
    Vector_y=cross(Vector_z, Vector_x);

    eje_referencia{x,1}=[Vector_x/norm(Vector_x); Vector_y/norm(Vector_y);Vector_z/norm(Vector_z)];
    
end

referencia=eje_referencia;
  
end