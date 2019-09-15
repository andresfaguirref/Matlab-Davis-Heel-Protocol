function centro_articular=Right_Knee_RC(names, dinamica, distancia)

X=Find_name(names,'r knee 1.X');
Y=Find_name(names,'r knee 1.Y');
Z=Find_name(names,'r knee 1.Z');

R_Knee1=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'r thigh.X');
Y=Find_name(names,'r thigh.Y');
Z=Find_name(names,'r thigh.Z');

R_Troch=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'r bar 1.X');
Y=Find_name(names,'r bar 1.Y');
Z=Find_name(names,'r bar 1.Z');

R_Bar1=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

punto_RC=[];

for x=1:length(R_Knee1(:,1))
    
    eje_ref = Eje_referencia(R_Knee1(x,:),R_Troch(x,:),R_Bar1(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),R_Knee1(x,:), -1*distancia/2); 
    
end

centro_articular = punto_RC;

end