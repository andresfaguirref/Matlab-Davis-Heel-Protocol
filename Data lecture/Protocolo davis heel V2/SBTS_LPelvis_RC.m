function centro_articular=SBTS_LPelvis_RC(names, rasis2, lasis2, static, punto)

%Me toca hallar el punto medio entre las espinas ilíacas con la distancia y
%mover punto

R_Asis=rasis2;

L_Asis=lasis2;

X=Find_name(names,'sacrum.X');
Y=Find_name(names,'sacrum.Y');
Z=Find_name(names,'sacrum.Z');

Sacrum=[static(:,X), static(:,Y), static(:,Z)];

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
    
punto_mueve_x(x,:)=Mover_punto(eje_ref(1,:),punto_nuevo(x,:), punto(1)); 
punto_mueve_y(x,:)=Mover_punto(eje_ref(2,:),punto_mueve_x(x,:), punto(2));
punto_mueve_z(x,:)=Mover_punto(eje_ref(3,:),punto_mueve_y(x,:), punto(3));
    
end

centro_articular = punto_mueve_z;

end