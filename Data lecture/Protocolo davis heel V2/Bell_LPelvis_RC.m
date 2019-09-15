function centro_articular=Bell_LPelvis_RC(names,rasis, lasis, dinamic)

%Me toca hallar el punto medio entre las espinas ilíacas con la distancia y
%mover punto


R_Asis=rasis;

L_Asis=lasis;

X=Find_name(names,'sacrum.X');
Y=Find_name(names,'sacrum.Y');
Z=Find_name(names,'sacrum.Z');

Sacrum=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

punto_nuevo=[];
vector=L_Asis-R_Asis;


for x=1:length(L_Asis(:,1))  
 
    distancia=norm(vector(x,:));
   
    punto_nuevo(x,:)=Mover_punto(vector(x,:),R_Asis(x,:), distancia/2); 
    
    
end

punto_mueve_x=[]; 
punto_mueve_y=[];
punto_mueve_z=[];

for x=1:length(L_Asis(:,1))

eje_ref = Eje_referencia(Sacrum(x,:),punto_nuevo(x,:),L_Asis(x,:));

distancia=norm(vector(x,:));
    
punto_mueve_x(x,:)=Mover_punto(eje_ref(1,:),punto_nuevo(x,:), -0.19*distancia); 
punto_mueve_y(x,:)=Mover_punto(eje_ref(2,:),punto_mueve_x(x,:), 0.36*distancia);
punto_mueve_z(x,:)=Mover_punto(eje_ref(3,:),punto_mueve_y(x,:), -0.3*distancia);
    
end

centro_articular = punto_mueve_z;

end