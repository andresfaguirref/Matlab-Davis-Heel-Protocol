function [referencia, Sternum]=SRefS_Trunk(names, static)

%Hallar el sistema de referencia del tronco

X=Find_name(names,'r should.X');
Y=Find_name(names,'r should.Y');
Z=Find_name(names,'r should.Z');

R_Should=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'l should.X');
Y=Find_name(names,'l should.Y');
Z=Find_name(names,'l should.Z');

L_Should=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'c7.X');
Y=Find_name(names,'c7.Y');
Z=Find_name(names,'c7.Z');

C7=[static(:,X), static(:,Y), static(:,Z)];

punto_nuevo=[];
vector=L_Should-R_Should;


for x=1:length(L_Should(:,1))  
 
    distancia=norm(vector(x,:));
   
    punto_nuevo(x,:)=Mover_punto(vector(x,:),R_Should(x,:), distancia/2); 
    
    
end

Sternum=punto_nuevo;
eje_referencia={};

for x=1:length(L_Should(:,1))
    
    Vector_z=L_Should(x,:)-punto_nuevo(x,:);
    Vector_int=C7(x,:)-punto_nuevo(x,:);

    Vector_x=cross(Vector_z,Vector_int);
    Vector_y=cross(Vector_z, Vector_x);

    eje_referencia{x,1}=[Vector_x/norm(Vector_x); Vector_y/norm(Vector_y);Vector_z/norm(Vector_z)];
    
end

referencia=eje_referencia;
  
end