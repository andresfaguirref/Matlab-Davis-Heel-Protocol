function [R_Asis2,L_Asis2]=Asis_correction(names, static,asis_distance)


%Utilizando el m�todo del 25% de distancia entre los trocanter 
%utilizando el eje de referencia de la rodilla

X=Find_name(names,'r asis.X');
Y=Find_name(names,'r asis.Y');
Z=Find_name(names,'r asis.Z');

R_Asis=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,'l asis.X');
Y=Find_name(names,'l asis.Y');
Z=Find_name(names,'l asis.Z');

L_Asis=[static(:,X), static(:,Y), static(:,Z)];
vector=L_Asis-R_Asis;
vector1=R_Asis-L_Asis;
punto_medio=[];

for x=1:length(L_Asis(:,1))  
 
    distancia=norm(vector(x,:));
    
    vector_unitario(x,:)=vector(x,:)/distancia;
    vector_unitario2(x,:)=vector1(x,:)/distancia;
    
    punto_medio(x,:)=Mover_punto(vector(x,:),R_Asis(x,:), distancia/2);  
       
end

Asis_dis=asis_distance/2;

for x=1:length(L_Asis(:,1))  
 
    nuevo_rasis(x,:)=Mover_punto(vector_unitario2(x,:),punto_medio(x,:), Asis_dis);
    nuevo_lasis(x,:)=Mover_punto(vector_unitario(x,:),punto_medio(x,:), Asis_dis);
       
end

R_Asis2=nuevo_rasis;
L_Asis2=nuevo_lasis;

end