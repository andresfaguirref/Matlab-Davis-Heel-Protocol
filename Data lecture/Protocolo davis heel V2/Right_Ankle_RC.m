function centro_articular=Right_Ankle_RC(names, dinamica, distancia)

X=Find_name(names,'r mall.X');
Y=Find_name(names,'r mall.Y');
Z=Find_name(names,'r mall.Z');

R_Mall=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'r knee 2.X');
Y=Find_name(names,'r knee 2.Y');
Z=Find_name(names,'r knee 2.Z');

R_Knee2=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'r bar 2.X');
Y=Find_name(names,'r bar 2.Y');
Z=Find_name(names,'r bar 2.Z');

R_Bar2=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

punto_RC=[];

for x=1:length(R_Knee2(:,1))
    
    eje_ref = Eje_referencia(R_Mall(x,:),R_Knee2(x,:),R_Bar2(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),R_Mall(x,:), -1*distancia/2); 
    
end

centro_articular = punto_RC;

end