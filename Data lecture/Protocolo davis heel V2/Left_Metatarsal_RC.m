function centro_articular=Left_Metatarsal_RC(names, dinamica, distancia)

X=Find_name(names,'l mall.X');
Y=Find_name(names,'l mall.Y');
Z=Find_name(names,'l mall.Z');

L_Mall=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l knee 2.X');
Y=Find_name(names,'l knee 2.Y');
Z=Find_name(names,'l knee 2.Z');

L_Knee2=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l bar 2.X');
Y=Find_name(names,'l bar 2.Y');
Z=Find_name(names,'l bar 2.Z');

L_Bar2=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l met.X');
Y=Find_name(names,'l met.Y');
Z=Find_name(names,'l met.Z');

L_5Met=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

punto_RC=[];

for x=1:length(L_Knee2(:,1))
    
    eje_ref = Eje_referencia(L_Mall(x,:),L_Knee2(x,:),L_Bar2(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),L_5Met(x,:), -1*distancia/2); 
    
end

centro_articular = punto_RC;

end