function [referencia, Sternum]=RefS_Trunk(names, dinamica)

%Hallar el sistema de referencia del tronco

X=Find_name(names,'r should.X');
Y=Find_name(names,'r should.Y');
Z=Find_name(names,'r should.Z');

R_Should=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l should.X');
Y=Find_name(names,'l should.Y');
Z=Find_name(names,'l should.Z');

L_Should=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'c7.X');
Y=Find_name(names,'c7.Y');
Z=Find_name(names,'c7.Z');

C7=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

punto_nuevo=[];
vector=L_Should-R_Should;


for x=1:length(L_Should(:,1))  
 
    distancia=norm(vector(x,:));
   
    punto_nuevo(x,:)=Mover_punto(vector(x,:),R_Should(x,:), distancia/2); 
    
    
end

Sternum=punto_nuevo;
eje_referencia={};

for x=1:length(L_Should(:,1))
    
    eje_ref = Eje_referencia(R_Should(x,:),L_Should(x,:),C7(x,:));
    eje_referencia{x}=eje_ref;
    
end

referencia=eje_referencia;
  
end