function centro_articular=DavisH_RPelvis_RC(names, dinamic, rasis, lasis,length_leg)
%Me toca hallar el punto medio entre las espinas ilíacas con la distancia y
%mover punto

R_Asis=rasis;

L_Asis=lasis;

X=Find_name(names,'sacrum.X');
Y=Find_name(names,'sacrum.Y');
Z=Find_name(names,'sacrum.Z');

Sacrum=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

punto_nuevo=[];
vector=L_Asis-R_Asis;


for x=1:length(L_Asis(:,1))  
 
    distancia=norm(vector(x,:));
   
    punto_nuevo(x,:)=Mover_punto(vector(x,:),R_Asis(x,:), distancia/2);     
    
end

%Calculamos las constantes que se van a utilizar para hallar las gráficas

theta=28.4*pi/180;
beta=18*pi/180;
xdis=0.1288*length_leg-0.04856;
radio=0.015/2;
C=0.115*length_leg-0.0153;

hip_y=(C*sin(theta)-0.5*norm(vector(1,:)));
hip_x=(-1*xdis-radio)*cos(beta)+C*cos(theta)*sin(beta);
hip_z=(-1*xdis-radio)*sin(beta)-C*cos(theta)*cos(beta);


punto_mueve_x=[]; 
punto_mueve_y=[];
punto_mueve_z=[];

for x=1:length(L_Asis(:,1))

eje_ref = Eje_referencia(Sacrum(x,:),punto_nuevo(x,:),L_Asis(x,:));

distancia=norm(vector(x,:));
    
punto_mueve_x(x,:)=Mover_punto(eje_ref(1,:),punto_nuevo(x,:), hip_x); 
punto_mueve_y(x,:)=Mover_punto(eje_ref(2,:),punto_mueve_x(x,:), hip_y);
punto_mueve_z(x,:)=Mover_punto(eje_ref(3,:),punto_mueve_y(x,:), hip_z);
    
end

centro_articular = punto_mueve_z;

end