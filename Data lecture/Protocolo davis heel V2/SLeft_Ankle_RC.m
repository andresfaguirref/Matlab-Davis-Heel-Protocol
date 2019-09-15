function centro_articular=SLeft_Ankle_RC(names, static, distancia)

X=Find_name(names,'l mall.X');
Y=Find_name(names,'l mall.Y');
Z=Find_name(names,'l mall.Z');

L_Mall=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'l knee 2.X');
Y=Find_name(names,'l knee 2.Y');
Z=Find_name(names,'l knee 2.Z');

L_Knee2=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'l bar 2.X');
Y=Find_name(names,'l bar 2.Y');
Z=Find_name(names,'l bar 2.Z');

L_Bar2=[static(:,X), static(:,Y), static(:,Z)];

punto_RC=[];

for x=1:length(L_Knee2(:,1))
    
    eje_ref = Eje_referencia(L_Mall(x,:),L_Knee2(x,:),L_Bar2(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),L_Mall(x,:), -1*distancia/2); 
    
end

centro_articular = punto_RC;

end