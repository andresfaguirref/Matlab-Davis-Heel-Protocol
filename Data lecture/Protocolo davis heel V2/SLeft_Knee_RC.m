function centro_articular=SLeft_Knee_RC(names, static, distancia)

X=Find_name(names,'l knee 1.X');
Y=Find_name(names,'l knee 1.Y');
Z=Find_name(names,'l knee 1.Z');

L_Knee1=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'l thigh.X');
Y=Find_name(names,'l thigh.Y');
Z=Find_name(names,'l thigh.Z');

L_Troch=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'l bar 1.X');
Y=Find_name(names,'l bar 1.Y');
Z=Find_name(names,'l bar 1.Z');

L_Bar1=[static(:,X), static(:,Y), static(:,Z)];

punto_RC=[];

for x=1:length(L_Knee1(:,1))
    
    eje_ref = Eje_referencia(L_Knee1(x,:),L_Troch(x,:),L_Bar1(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),L_Knee1(x,:), -1*distancia/2); 
    
end

centro_articular = punto_RC;

end