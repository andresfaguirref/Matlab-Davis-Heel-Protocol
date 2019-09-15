function centro_articular=PerC25_RPelvis_RC(names, dinamica)

%teniendo en cuenta el vector formado entre los dos trocanter en cada
%instante de tiempo, y el 25% de esta distancia para ubicar el nuevo punto.

X=Find_name(names,'r thigh.X');
Y=Find_name(names,'r thigh.Y');
Z=Find_name(names,'r thigh.Z');

R_Troch=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];

X=Find_name(names,'l thigh.X');
Y=Find_name(names,'l thigh.Y');
Z=Find_name(names,'l thigh.Z');

L_Troch=[dinamica(:,X), dinamica(:,Y), dinamica(:,Z)];


punto_nuevo=[];
dis_troch=L_Troch-R_Troch;

punto_RC=[];

for x=1:length(L_Troch(:,1))

distancia=norm(dis_troch(x,:));

vector_troch=(L_Troch(x,:)-R_Troch(x,:))/distancia;
    
punto_RC(x,:)=Mover_punto(vector_troch,R_Troch(x,:), 0.25*distancia); 
    
end

centro_articular = punto_RC;

end