function centro_articular=SRight_Metatarsal_RC(names, static, distancia)

X=Find_name(names,'r mall.X');
Y=Find_name(names,'r mall.Y');
Z=Find_name(names,'r mall.Z');

R_Mall=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'r knee 2.X');
Y=Find_name(names,'r knee 2.Y');
Z=Find_name(names,'r knee 2.Z');

R_Knee2=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'r bar 2.X');
Y=Find_name(names,'r bar 2.Y');
Z=Find_name(names,'r bar 2.Z');

R_Bar2=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'r met.X');
Y=Find_name(names,'r met.Y');
Z=Find_name(names,'r met.Z');

R_5Met=[static(:,X), static(:,Y), static(:,Z)];

punto_RC=[];

for x=1:length(R_Knee2(:,1))
    
    eje_ref = Eje_referencia(R_Mall(x,:),R_Knee2(x,:),R_Bar2(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),R_5Met(x,:), -1*distancia/2); 
    
end

centro_articular = punto_RC;

end