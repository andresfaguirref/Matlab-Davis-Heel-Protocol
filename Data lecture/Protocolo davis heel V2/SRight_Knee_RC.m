function centro_articular=SRight_Knee_RC(names, static, distancia)

X=Find_name(names,'r knee 1.X');
Y=Find_name(names,'r knee 1.Y');
Z=Find_name(names,'r knee 1.Z');

R_Knee1=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'r thigh.X');
Y=Find_name(names,'r thigh.Y');
Z=Find_name(names,'r thigh.Z');

R_Troch=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'r bar 1.X');
Y=Find_name(names,'r bar 1.Y');
Z=Find_name(names,'r bar 1.Z');

R_Bar1=[static(:,X), static(:,Y), static(:,Z)];

punto_RC=[];

for x=1:length(R_Knee1(:,1))
    
    eje_ref = Eje_referencia(R_Knee1(x,:),R_Troch(x,:),R_Bar1(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),R_Knee1(x,:), -1*distancia/2); 
    
end

centro_articular = punto_RC;

end