function centro_articular=PC25_LPelvis_RC(names, dinamica)

%Utilizando el método del 25% de distancia entre los trocanter 
%utilizando el eje de referencia de la rodilla

X=Find_name(names,'l knee 1.X');
Y=Find_name(names,'l knee 1.Y');
Z=Find_name(names,'l knee 1.Z');

L_Knee1=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l thigh.X');
Y=Find_name(names,'l thigh.Y');
Z=Find_name(names,'l thigh.Z');

L_Troch=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l bar 1.X');
Y=Find_name(names,'l bar 1.Y');
Z=Find_name(names,'l bar 1.Z');

L_Bar1=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'r thigh.X');
Y=Find_name(names,'r thigh.Y');
Z=Find_name(names,'r thigh.Z');

R_Troch=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l thigh.X');
Y=Find_name(names,'l thigh.Y');
Z=Find_name(names,'l thigh.Z');

L_Troch=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

dis_troch=L_Troch-R_Troch;

punto_RC=[];

for x=1:length(L_Knee1(:,1))
    
    eje_ref = Eje_referencia(L_Knee1(x,:),L_Troch(x,:),L_Bar1(x,:));

    distancia=norm(dis_troch(x,:));
    
    punto_RC(x,:)=Mover_punto(eje_ref(2,:),L_Troch(x,:), -0.25*distancia); 
    
end

centro_articular = punto_RC;

end